#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/MathFunctions.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/GeometryFunctions.h>
#include <raytracer/geometry/Ray.h>
#include <raytracer/geometry/MeshFunction.h>
#include <raytracer/physics/LaserPropagation.h>
#include <raytracer/physics/GradientCalculators.h>
#include <chrono>
#include "mfem.hpp"

using namespace raytracer::geometry;
using namespace raytracer::physics;

double density(const mfem::Vector &x) {
    return 12.8e20 * x(0);
}

struct EndAbsorber {
    explicit EndAbsorber(MeshFunction &absorbedEnergy) : absorbedEnergy(absorbedEnergy) {}

    void operator()(const LaserRay &laserRay) {
        const auto &intersection = laserRay.intersections.back();
        const auto element = intersection.previousElement;
        if (!element) return;
        absorbedEnergy.addValue(*element, laserRay.energy.asDouble);
    }

private:
    MeshFunction &absorbedEnergy;
};





int main(int, char *[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    //DiscreteLine side{};
    //side.segmentCount = 100;
    //side.length = 1;
    //auto mfemMesh = constructRectangleMesh(side, side);

    auto mfemMesh = std::make_unique<mfem::Mesh>("test_mesh.vtk", 1, 0);


    mfem::L2_FECollection l2FiniteElementCollection(0, 2);
    mfem::H1_FECollection h1FiniteElementCollection(1, 2);
    mfem::FiniteElementSpace l2FiniteElementSpace(mfemMesh.get(), &l2FiniteElementCollection);
    mfem::FiniteElementSpace h1FiniteElementSpace(mfemMesh.get(), &h1FiniteElementCollection);

    mfem::GridFunction densityGridFunction(&l2FiniteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    mfem::GridFunction absorbedEnergyGridFunction(&l2FiniteElementSpace);
    absorbedEnergyGridFunction = 0;
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());
    MfemMeshFunction densityMeshFunction(densityGridFunction, l2FiniteElementSpace);
    MfemMeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, l2FiniteElementSpace);

    //ConstantGradientCalculator gradientCalculator(Vector(12.8e20, 0));
    H1GradientCalculator gradientCalculator(l2FiniteElementSpace, h1FiniteElementSpace);
    gradientCalculator.updateDensity(densityGridFunction);

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point) { return Vector(1, 0.3); },
            Gaussian(0.1),
            Point(-1.1, -0.4),
            Point(-1.1, -0.8)
    );

    laser.generateRays(100);
    laser.generateIntersections(
            mesh, SnellsLaw(densityMeshFunction, gradientCalculator),
            DontStop());

    EndAbsorber endAbsorber(absorbedEnergyMeshFunction);
    for (const auto &laserRay : laser.getRays()) {
        endAbsorber(laserRay);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto nanosecondDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    auto duration = nanosecondDuration * 1e-9;

    std::cout << "Execution took: " << duration << "s." << std::endl;

    laser.saveRaysToJson("rays.json");

    //GLVIS-------------------------------------------------------------------------------------------------------------
    char vishost[] = "localhost";
    int visport = 19916;
    mfem::socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    //sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
    sol_sock << "solution\n" << *mfemMesh << densityGridFunction << std::flush;
}