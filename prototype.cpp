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

double densityFunction(const mfem::Vector &x) {
    return 12.8e26 * x(0);
}

double temperatureFunction(const mfem::Vector &) {
    return 20;
}

double ionizationFunction(const mfem::Vector &) {
    return 0;
}

class StepGradient : public GradientCalculator {
public:
    Vector getGradient(const Intersection &intersection) const override {
        return intersection.face->getNormal();
    }
};

class AbsorptionModel {
public:
    virtual Energy getEnergyChange(
            const Intersection &previousIntersection,
            const Intersection &currentIntersection,
            const Energy &currentEnergy,
            const LaserRay& laserRay
    ) const = 0;
};

class EnergyAbsorber {
public:
    void addModel(const AbsorptionModel* model) {
        models.emplace_back(model);
    }

    void absorb(const Laser &laser, MeshFunction &absorbedEnergy) {
        for (const auto &laserRay : laser.getRays()) {
            this->absorbLaserRay(laserRay, absorbedEnergy);
        }
    }

private:
    std::vector<const AbsorptionModel*> models{};

    void absorbLaserRay(const LaserRay &laserRay, MeshFunction &absorbedEnergy) {
        const auto &intersections = laserRay.intersections;
        auto intersectionIt = std::next(std::begin(intersections));
        auto previousIntersectionIt = std::begin(intersections);

        auto currentEnergy = laserRay.energy.asDouble;

        for (; intersectionIt != std::end(intersections); ++intersectionIt, ++previousIntersectionIt) {
            for (const auto &model : this->models) {
                auto absorbed = model->getEnergyChange(
                        *previousIntersectionIt,
                        *intersectionIt,
                        Energy{currentEnergy},
                        laserRay
                ).asDouble;
                currentEnergy -= absorbed;
                absorbedEnergy.addValue(*(intersectionIt->previousElement), absorbed);
            }
        }
    }
};

struct BremsstrahlungModel : public AbsorptionModel {
    explicit BremsstrahlungModel(
            const MeshFunction &density,
            const MeshFunction &temperature,
            const MeshFunction &ionization
    ) :
            _density(density),
            _temperature(temperature),
            _ionization(ionization) {}

    Energy getEnergyChange(
            const Intersection &previousIntersection,
            const Intersection &currentIntersection,
            const Energy &currentEnergy,
            const LaserRay& laserRay
            ) const override {
        const auto &element = currentIntersection.previousElement;
        if (!element) return Energy{0};
        const auto &previousPoint = previousIntersection.orientation.point;
        const auto &point = currentIntersection.orientation.point;

        const auto distance = (point - previousPoint).getNorm();
        const auto density = Density{this->_density.getValue(*element)};
        const auto temperature = Temperature{this->_temperature.getValue(*element)};
        const auto ionization = this->_ionization.getValue(*element);

        SpitzerFrequencyCalculator frequencyCalculator;
        auto frequency = frequencyCalculator.getCollisionalFrequency(density, temperature, laserRay.wavelength,
                                                                     ionization);
        const auto exponent = -laserRay.getInverseBremsstrahlungCoeff(density, frequency) * distance;

        auto newEnergy = currentEnergy.asDouble * std::exp(exponent);
        return Energy{currentEnergy.asDouble - newEnergy};
    }

private:
    const MeshFunction &_density;
    const MeshFunction &_temperature;
    const MeshFunction &_ionization;
};

int main(int, char *[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    //DiscreteLine side{};
    //side.segmentCount = 100;
    //side.length = 1e-6;
    //auto mfemMesh = constructRectangleMesh(side, side);

    auto mfemMesh = std::make_unique<mfem::Mesh>("test_mesh_small.vtk", 1, 0);


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

    //ConstantGradientCalculator gradientCalculator(Vector(12.8e20, 0));
    //H1GradientCalculator gradientCalculator(l2FiniteElementSpace, h1FiniteElementSpace);
    StepGradient gradientCalculator;
    //gradientCalculator.updateDensity(densityGridFunction);

    SpitzerFrequencyCalculator spitzerFrequencyCalculator;
    SnellsLaw snellsLaw(
            densityMeshFunction,
            temperatureMeshFunction,
            ionizationMeshFunction,
            gradientCalculator,
            spitzerFrequencyCalculator
    );

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point) { return Vector(1, 0.3); },
            Gaussian(0.1),
            Point(-0.1e-6, 0.5e-6),
            Point(-0.1e-6, 0)
    );

    laser.generateRays(100);
    laser.generateIntersections(
            mesh, snellsLaw,
            DontStop());

    EnergyAbsorber absorber;
    BremsstrahlungModel bremsstrahlungModel(densityMeshFunction, temperatureMeshFunction, ionizationMeshFunction);
    absorber.addModel(&bremsstrahlungModel);
    absorber.absorb(laser, absorbedEnergyMeshFunction);

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
    sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
    //sol_sock << "solution\n" << *mfemMesh << densityGridFunction << std::flush;
}