#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/MathFunctions.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/GeometryFunctions.h>
#include <raytracer/geometry/Ray.h>
#include <raytracer/geometry/MeshFunction.h>
#include <chrono>
#include "mfem.hpp"

using namespace raytracer::geometry;
using namespace raytracer::physics;

struct StopAtCritical {
    StopAtCritical(const LaserRay &laserRay, const MeshFunction &density) :
            laserRay(laserRay), density(density) {}

    bool operator()(const Intersection &intersection) {
        const auto element = *intersection.element;

        auto currentDensity = density[element];
        auto criticalDensity = laserRay.getCriticalDensity();
        return currentDensity > criticalDensity.asDouble;
    }

private:
    const LaserRay &laserRay;
    const MeshFunction &density;
};


std::unique_ptr<Intersection> continueStraight(const Intersection &intersection) {
    const auto element = intersection.element;
    return findClosestIntersection(intersection.orientation, element->getFaces());
}


double density(const mfem::Vector &x) {
    return 12.8e20 * x(0);
}

void absorbEnergyAtRayEnd(const LaserRay &laserRay, MeshFunction& absorbedEnergy) {
    const auto &intersection = laserRay.intersections.back();
    const auto element = *intersection.element;
    absorbedEnergy[element] += laserRay.energy.asDouble;
}

int main(int argc, char *argv[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    bool visualization = true;

    DiscreteLine side{};
    side.segmentCount = 100;
    side.width = 1;
    auto mfemMesh = constructRectangleMesh(side, side);

    mfem::L2_FECollection finiteElementCollection(0, 2);
    mfem::FiniteElementSpace finiteElementSpace(mfemMesh.get(), &finiteElementCollection);

    mfem::GridFunction densityGridFunction(&finiteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    mfem::GridFunction absorbedEnergyGridFunction(&finiteElementSpace);
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());
    MeshFunction densityMeshFunction(densityGridFunction, finiteElementSpace);
    MeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, finiteElementSpace);

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point &point) { return Vector(1, 0.1); },
            Gaussian(0.2),
            Point(-1, 0.3),
            Point(-1, 0.7)
    );

    auto laserRays = laser.generateRays(10000);

    for (auto &laserRay : laserRays) {
        laserRay.generateIntersections(mesh, continueStraight, StopAtCritical(laserRay, densityMeshFunction));
    }

    for (const auto &laserRay : laserRays) {
        absorbEnergyAtRayEnd(laserRay, absorbedEnergyMeshFunction);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto nanosecondDuration = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
    auto duration = nanosecondDuration * 1e-9;

    std::cout << "Execution took: " << duration << "s." << std::endl;


    //GLVIS-------------------------------------------------------------------------------------------------------------
    char vishost[] = "localhost";
    int visport = 19916;
    mfem::socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
}