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
    explicit StopAtCritical(const MeshFunction &density) :
    density(density) {}

    bool operator()(const Intersection &intersection, const LaserRay &laserRay) {
        const auto element = *intersection.nextElement;

        auto currentDensity = density[element];
        auto criticalDensity = laserRay.getCriticalDensity();
        return currentDensity > criticalDensity.asDouble;
    }

private:
    const MeshFunction &density;
};


std::unique_ptr<Intersection> continueStraight(const Intersection &intersection) {
    const auto element = intersection.nextElement;
    auto found = findClosestIntersection(intersection.orientation, element->getFaces(), intersection.face);
    return found;
}

struct SnellsLaw {
    explicit SnellsLaw(const MeshFunction &density) :
            density(density) {}

    std::unique_ptr<Intersection> operator()(const Intersection &intersection) {
        const auto previousElement = intersection.previousElement;
        const auto nextElement = intersection.nextElement;

        if (!previousElement) {
            return findClosestIntersection(intersection.orientation, nextElement->getFaces(), intersection.face);
        }

        const double n1 = std::sqrt(1 - density[*previousElement] / 6e20);
        const double n2 = std::sqrt(1 - density[*nextElement] / 6e20);

        const Vector gradient(12.8e20, 0);
        const auto& direction = intersection.orientation.direction;
        const auto l = 1/direction.getNorm() * direction;
        auto n = 1/gradient.getNorm() * gradient;
        auto c = (-1) * n * l;
        if (c < 0){
            c = -c;
            n = (-1) * n;
        }
        const double r = n1 / n2;

        Vector newDirection{};
        auto root = 1 - r*r * (1 - c*c);
        if (root > 0) {
            newDirection = r * l + (r*c - std::sqrt(root))*n;
        } else {
            newDirection = l + 2*c*n;
        }

        auto newIntersection = findClosestIntersection(
                intersection.orientation,
                nextElement->getFaces(),
                intersection.face);
        newIntersection->orientation.direction = newDirection;
        return newIntersection;
    }
private:
    const MeshFunction &density;
};


double density(const mfem::Vector &x) {
    return 12.8e20 * x(0);
}

struct EndAbsorber {
    explicit EndAbsorber(MeshFunction& absorbedEnergy): absorbedEnergy(absorbedEnergy) {}

    void operator()(const LaserRay &laserRay) {
        const auto &intersection = laserRay.intersections.back();
        const auto element = intersection.nextElement;
        if (!element) return;
        absorbedEnergy[*element] += laserRay.energy.asDouble;
    }
private:
    MeshFunction& absorbedEnergy;
};

int main(int argc, char *argv[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    DiscreteLine side{};
    side.segmentCount = 100;
    side.length = 1;
    //auto mfemMesh = constructRectangleMesh(side, side);

    auto mfemMesh = std::make_unique<mfem::Mesh>("test_mesh.vtk", 1, 0);


    mfem::L2_FECollection finiteElementCollection(0, 2);
    mfem::FiniteElementSpace finiteElementSpace(mfemMesh.get(), &finiteElementCollection);

    mfem::GridFunction densityGridFunction(&finiteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    mfem::GridFunction absorbedEnergyGridFunction(&finiteElementSpace);
    absorbedEnergyGridFunction = 0;
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());
    MeshFunction densityMeshFunction(densityGridFunction, finiteElementSpace);
    MeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, finiteElementSpace);

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point &point) { return Vector(1, 0); },
            Gaussian(0.4),
            Point(-1.1, 0.8),
            Point(-1.1, -0.8)
    );

    laser.generateRays(1000);
    laser.generateIntersections(mesh, continueStraight, StopAtCritical(densityMeshFunction));

    EndAbsorber endAbsorber(absorbedEnergyMeshFunction);
    for (const auto &laserRay : laser.getRays()) {
        endAbsorber(laserRay);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto nanosecondDuration = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
    auto duration = nanosecondDuration * 1e-9;

    std::cout << "Execution took: " << duration << "s." << std::endl;

    laser.saveRaysToJson("rays.json");

    //GLVIS-------------------------------------------------------------------------------------------------------------
    char vishost[] = "localhost";
    int visport = 19916;
    mfem::socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
}