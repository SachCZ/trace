#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/MeshFunction.h>
#include <raytracer/physics/CollisionalFrequency.h>
#include <raytracer/physics/Gradient.h>
#include <raytracer/physics/Refraction.h>
#include <raytracer/physics/Propagation.h>
#include <raytracer/physics/Termination.h>
#include <raytracer/physics/Absorption.h>
#include "mfem.hpp"
#include "yaml-cpp/yaml.h"

raytracer::Vector parseVector(const YAML::Node& node){
    return {node["x"].as<double>(), node["y"].as<double>()};
}

raytracer::Point parsePoint(const YAML::Node& node){
    return {node["x"].as<double>(), node["y"].as<double>()};
}

int main(int, char * argv[]) {
    using namespace raytracer;

    YAML::Node config = YAML::LoadFile(argv[1]);

    std::string meshFilename = config["mesh_file"].as<std::string>();
    std::string densityFilename = config["density_file"].as<std::string>();
    std::string temperatureFilename = config["temperature_file"].as<std::string>();
    std::string ionizationFilename = config["ionization_file"].as<std::string>();
    std::string raysOutputFilename = "output/rays.msgpack";
    std::string absorbedEnergyFilename = "output/absorbed_energy.gf";

    auto mfemMesh = std::make_unique<mfem::Mesh>(meshFilename.c_str(), 1, 1, false);
    MfemMesh mesh(mfemMesh.get());

    mfem::L2_FECollection l2FiniteElementCollection(0, 2);
    mfem::H1_FECollection h1FiniteElementCollection(1, 2);
    mfem::FiniteElementSpace l2FiniteElementSpace(mfemMesh.get(), &l2FiniteElementCollection);
    mfem::FiniteElementSpace h1FiniteElementSpace(mfemMesh.get(), &h1FiniteElementCollection);

    mfem::GridFunction absorbedEnergyGridFunction(&l2FiniteElementSpace);
    absorbedEnergyGridFunction = 0;
    MfemMeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, l2FiniteElementSpace);

    std::ifstream temperatureFile(temperatureFilename);
    mfem::GridFunction temperatureGridFunction(mfemMesh.get(), temperatureFile);
    MfemMeshFunction temperatureMeshFunction(temperatureGridFunction, l2FiniteElementSpace);

    std::ifstream ionizationFile(ionizationFilename);
    mfem::GridFunction ionizationGridFunction(mfemMesh.get(), ionizationFile);
    MfemMeshFunction ionizationMeshFunction(ionizationGridFunction, l2FiniteElementSpace);

    std::ifstream densityFile(densityFilename);
    mfem::GridFunction densityGridFunction(mfemMesh.get(), densityFile);
    double massUnit = 1.6605e-24;
    double A = 55.845;
    mfem::Vector electronDensity(densityGridFunction.Size());
    densityGridFunction.GetTrueDofs(electronDensity);
    mfem::Vector ionization(ionizationGridFunction.Size());
    ionizationGridFunction.GetTrueDofs(ionization);
    for (int i = 0; i < densityGridFunction.Size(); i++){
        electronDensity[i] = electronDensity[i] * ionization[i] / massUnit / A;
    }
    densityGridFunction.SetFromTrueDofs(electronDensity);

    MfemMeshFunction electronDensityMeshFunction(densityGridFunction, l2FiniteElementSpace);

    IntersectionSet allIntersections;

    for (const auto& node : config["lasers"]){
        Length wavelength{node["wavelength"].as<double>()};
        auto direction{parseVector(node["direction"])};
        auto spatialFWHM = node["spatial_FWHM"].as<double>();
        Energy energy{node["energy"].as<double>()};
        Point startPoint(parsePoint(node["start_point"]));
        Point endPoint(parsePoint(node["end_point"]));
        auto raysCount = node["rays_count"].as<int>();
        auto estimateBremsstrahlung = node["estimate_bremsstrahlung"].as<bool>();
        auto estimateResonance = node["estimate_resonance"].as<bool>();
        auto estimateGain = node["estimate_gain"].as<bool>();
        auto gradientType = node["gradient_type"].as<std::string>();

        Laser laser(
                wavelength,
                [&direction](const Point &) { return direction; },
                Gaussian(spatialFWHM, energy.asDouble, 0),
                startPoint,
                endPoint,
                raysCount
        );

        std::unique_ptr<Gradient> gradient;
        auto h1Function = projectL2toH1(densityGridFunction, l2FiniteElementSpace, h1FiniteElementSpace);

        if (gradientType == "ls"){
            gradient = std::make_unique<LinearInterpolation>(getHouseholderGradientAtPoints(mesh, electronDensityMeshFunction));
        } else if (gradientType == "h1") {
            gradient = std::make_unique<H1Gradient>(h1Function, *mfemMesh);
        }

        SpitzerFrequency spitzerFrequency;
        ColdPlasma coldPlasma;
        Marker reflected;
        SnellsLaw snellsLaw(
                electronDensityMeshFunction,
                temperatureMeshFunction,
                ionizationMeshFunction,
                *gradient,
                spitzerFrequency,
                coldPlasma,
                laser.wavelength,
                &reflected
        );


        auto initialDirections = generateInitialDirections(laser);
        auto intersections = generateIntersections(
                mesh,
                initialDirections,
                snellsLaw,
                intersectStraight,
                DontStop()
        );
        allIntersections.insert(allIntersections.end(), intersections.begin(), intersections.end());

        AbsorptionController absorber;

        std::unique_ptr<Bremsstrahlung> bremsstrahlungModel;
        if (estimateBremsstrahlung){
            bremsstrahlungModel = std::make_unique<Bremsstrahlung>(
                    electronDensityMeshFunction,
                    temperatureMeshFunction,
                    ionizationMeshFunction,
                    spitzerFrequency,
                    coldPlasma,
                    laser.wavelength
            );
            absorber.addModel(bremsstrahlungModel.get());
        }

        std::unique_ptr<Resonance> resonance;
        ClassicCriticalDensity classicCriticalDensity;
        if (estimateResonance){
            resonance = std::make_unique<Resonance>(*gradient, classicCriticalDensity, laser.wavelength, reflected);
            absorber.addModel(resonance.get());
        }

        std::unique_ptr<XRayGain> gain;
        std::unique_ptr<MfemMeshFunction> gainMeshFunction;
        std::unique_ptr<mfem::GridFunction> gainGridFunction;
        if (estimateGain){
            std::ifstream gainFile(config["gain_file"].as<std::string>());
            gainGridFunction = std::make_unique<mfem::GridFunction>(mfemMesh.get(), gainFile);
            gainMeshFunction = std::make_unique<MfemMeshFunction>(*gainGridFunction, l2FiniteElementSpace);
            gain = std::make_unique<XRayGain>(*gainMeshFunction);
            absorber.addModel(gain.get());
        }

        std::cout << stringifyAbsorptionSummary(
                absorber.absorb(intersections, generateInitialEnergies(laser), absorbedEnergyMeshFunction)
        );
    }

    std::ofstream raysFile(raysOutputFilename);
    raysFile << stringifyRaysToMsgpack(allIntersections);
    std::ofstream absorbedResult(absorbedEnergyFilename);
    absorbedResult << absorbedEnergyMeshFunction;
}