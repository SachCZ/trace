#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/MeshFunction.h>
#include <chrono>
#include <raytracer/physics/CollisionalFrequency.h>
#include <raytracer/physics/Gradient.h>
#include <raytracer/physics/Refraction.h>
#include <raytracer/physics/Propagation.h>
#include <raytracer/physics/Termination.h>
#include <raytracer/physics/Absorption.h>
#include <set>
#include "mfem.hpp"

using namespace raytracer;

double densityFunction(const mfem::Vector &x) {
    return 6.3e25 * x(0) + 12.8e20 / 2;
    //return 12.8e26 * x(0);
}

double temperatureFunction(const mfem::Vector &) {
    return 28;
}

double ionizationFunction(const mfem::Vector &) {
    return 22;
}

class Resonance : public AbsorptionModel {
public:
    Resonance(const Gradient &gradientCalculator, const Marker& reflectedMarker):
            gradientCalculator(gradientCalculator), reflectedMarker(reflectedMarker) {}

    Energy getEnergyChange(
            const Intersection &previousIntersection,
            const Intersection &currentIntersection,
            const Energy &currentEnergy,
            const LaserRay &laserRay
    ) const override {
        if (!Resonance::isResonating(*currentIntersection.previousElement)) return Energy{0};

        auto grad = gradientCalculator.get(
                currentIntersection.pointOnFace,
                *currentIntersection.previousElement,
                *currentIntersection.nextElement
        );
        auto dir = (currentIntersection.pointOnFace.point - previousIntersection.pointOnFace.point);
        auto q = Resonance::getQ(laserRay, dir, grad);
        auto term = q*std::exp(-4.0/3.0*std::pow(q, 3.0/2.0))/(q + 0.48) * M_PI / 2.0;
        return Energy{currentEnergy.asDouble*term};
    }

private:
    const Gradient &gradientCalculator;
    const Marker& reflectedMarker;

    bool isResonating(const Element& element) const {
        return reflectedMarker.isMarked(element);
    }

    static double getQ(const LaserRay &laserRay, Vector dir, Vector grad) {
        auto dir_norm = dir.getNorm();
        auto grad_norm = grad.getNorm();
        auto lamb = laserRay.wavelength.asDouble;
        auto ne_crit = laserRay.getCriticalDensity().asDouble;
        auto sin2phi = 1 - std::pow(grad * dir / grad_norm / dir_norm, 2);
        return std::pow(2 * M_PI / lamb * ne_crit / grad_norm, 2.0 / 3.0) * sin2phi;
    }
};

int main(int, char *[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    //DiscreteLine side{};
    //side.segmentCount = 100;
    //side.length = 1e-6;
    //auto mfemMesh = constructRectangleMesh(side, side);

    auto mfemMesh = std::make_unique<mfem::Mesh>("micrometr_mesh.vtk", 1, 0);


    mfem::L2_FECollection l2FiniteElementCollection(0, 2);
    mfem::H1_FECollection h1FiniteElementCollection(1, 2);
    mfem::FiniteElementSpace l2FiniteElementSpace(mfemMesh.get(), &l2FiniteElementCollection);
    mfem::FiniteElementSpace h1FiniteElementSpace(mfemMesh.get(), &h1FiniteElementCollection);

    mfem::GridFunction absorbedEnergyGridFunction(&l2FiniteElementSpace);
    absorbedEnergyGridFunction = 0;
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());

    //Density
    mfem::GridFunction densityGridFunction(&l2FiniteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(densityFunction);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);
    MfemMeshFunction densityMeshFunction(densityGridFunction, l2FiniteElementSpace);

    mfem::GridFunction temperatureGridFunction(&l2FiniteElementSpace);
    mfem::FunctionCoefficient temperatureFunctionCoefficient(temperatureFunction);
    temperatureGridFunction.ProjectCoefficient(temperatureFunctionCoefficient);
    MfemMeshFunction temperatureMeshFunction(temperatureGridFunction, l2FiniteElementSpace);

    mfem::GridFunction ionizationGridFunction(&l2FiniteElementSpace);
    mfem::FunctionCoefficient ionizationFunctionCoefficient(ionizationFunction);
    ionizationGridFunction.ProjectCoefficient(ionizationFunctionCoefficient);
    MfemMeshFunction ionizationMeshFunction(ionizationGridFunction, l2FiniteElementSpace);

    MfemMeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, l2FiniteElementSpace);

    //ConstantGradient constantGradient(Vector(6.3e25, 0));
    H1Gradient h1Gradient(l2FiniteElementSpace, h1FiniteElementSpace);
    //StepGradient stepGradient;
    h1Gradient.updateDensity(densityGridFunction);

    SpitzerFrequency spitzerFrequency;

    Marker reflectedMarker;
    SnellsLaw snellsLaw(
            densityMeshFunction,
            temperatureMeshFunction,
            ionizationMeshFunction,
            h1Gradient,
            spitzerFrequency,
            &reflectedMarker
    );

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point) { return Vector(1, 1); },
            Gaussian(0.3e-5), // [](double){return 1;}
            Point(-0.51e-5, -0.3e-5),
            Point(-0.51e-5, -0.5e-5)
    );

    laser.generateRays(100);

    laser.generateIntersections(mesh, snellsLaw, intersectStraight, DontStop());

    AbsorptionController absorber;
    Bremsstrahlung bremsstrahlungModel(
            densityMeshFunction,
            temperatureMeshFunction,
            ionizationMeshFunction,
            spitzerFrequency
    );
    Resonance resonance(h1Gradient, reflectedMarker);
    absorber.addModel(&bremsstrahlungModel);
    absorber.addModel(&resonance);
    absorber.absorb(laser, absorbedEnergyMeshFunction);

    auto end = std::chrono::high_resolution_clock::now();
    auto nanosecondDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    auto duration = nanosecondDuration * 1e-9;

    std::cout << "Execution took: " << duration << "s." << std::endl;

    laser.saveRaysToJson("rays.json");

    std::ofstream absorbedResult("absorbedEnergy.txt");
    absorbedEnergyGridFunction.Save(absorbedResult);

    //GLVIS-------------------------------------------------------------------------------------------------------------
    char vishost[] = "localhost";
    int visport = 19916;
    mfem::socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
    //sol_sock << "solution\n" << *mfemMesh << densityGridFunction << std::flush;
}