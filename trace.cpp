#include <raytracer.h>
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

    auto meshFilename = config["mesh_file"].as<std::string>();
    auto temperatureFilename = config["temperature_file"].as<std::string>();
    auto ionizationFilename = config["ionization_file"].as<std::string>();
    std::string raysOutputFilename = "output/rays.msgpack";
    std::string absorbedEnergyFilename = "output/absorbed_energy.gf";

    MfemMesh mesh(meshFilename);

    MfemL20Space space(mesh);
    //mfem::H1_FECollection h1FiniteElementCollection(1, 2);
    //mfem::FiniteElementSpace h1FiniteElementSpace(mfemMesh.get(), &h1FiniteElementCollection);

    MfemMeshFunction absorbedEnergy(space, [](Point){return 0;});

    std::ifstream temperatureFile(temperatureFilename);
    MfemMeshFunction temperature(space, temperatureFile);

    std::ifstream ionizationFile(ionizationFilename);
    MfemMeshFunction ionization(space, ionizationFile);

    std::unique_ptr<MfemMeshFunction> eleDens;
    if (config["density_file"]){
        std::ifstream densityFile(config["density_file"].as<std::string>());
        MfemMeshFunction density(space, densityFile);
        eleDens = std::make_unique<MfemMeshFunction>(space, [&](const Element& e) {
            double massUnit = 1.6605e-24;
            double A = 55.845;
            return density.getValue(e) * ionization.getValue(e) / massUnit / A;
        });
    } else if (config["ele_dens_profile"] && config["max_ele_dens"]){
        if (config["ele_dens_profile"].as<std::string>() == "parabolic"){
            eleDens = std::make_unique<MfemMeshFunction>(space, [&](const Point& point) {
                auto maxDens = config["max_ele_dens"].as<double>();
                return maxDens * (1 - std::pow(point.x, 2));
            });
        } else {
            throw std::logic_error("unknown profile");
        }
    } else {
        throw std::logic_error("Invalid density spec");
    }

    IntersectionSet allIntersections;

    for (const auto& node : config["lasers"]){
        auto gradientType = node["gradient_type"].as<std::string>();

        VectorField gradAtPoints;
        if (gradientType == "ls"){
            gradAtPoints = calcHousGrad(mesh, *eleDens);
        } else if (gradientType == "h1") {
            //gradient = std::make_unique<H1Gradient>(h1Function, *mfemMesh);
        }
        auto gradient = LinInterGrad(gradAtPoints);

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

        Laser laser{
                wavelength,
                [&direction](const Point &) { return direction; },
                Gaussian(spatialFWHM, energy.asDouble, 0),
                startPoint,
                endPoint,
                raysCount
        };

        MfemMeshFunction frequency (space, [&](const Element& e) {
            return calcSpitzerFreq(eleDens->getValue(e), temperature.getValue(e), ionization.getValue(e), wavelength);
        });
        MfemMeshFunction refractIndex(space, [&](const Element& e){
            return calcRefractIndex(eleDens->getValue(e), Length{wavelength}, frequency.getValue(e));
        });
        Marker reflected;
        SnellsLaw snellsLaw(gradient, refractIndex, &reflected);


        auto initialDirections = generateInitialDirections(laser);
        auto intersections = findIntersections(mesh, initialDirections, snellsLaw, intersectStraight, dontStop);
        allIntersections.insert(allIntersections.end(), intersections.begin(), intersections.end());

        EnergyExchangeController exchangeController;

        MfemMeshFunction invBremssCoeff(space, [&](const Element& e){
            return calcInvBremssCoeff(eleDens->getValue(e), wavelength, frequency.getValue(e));
        });

        Bremsstrahlung bremss(invBremssCoeff);
        if (estimateBremsstrahlung){
            exchangeController.addModel(&bremss);
        }

        Resonance resonance(gradient, laser.wavelength, reflected);
        if (estimateResonance){
            exchangeController.addModel(&resonance);
        }

        std::unique_ptr<MfemMeshFunction> gain;
        std::unique_ptr<XRayGain> xRayGain;
        if (estimateGain){
            std::ifstream gainFile(config["gain_file"].as<std::string>());
            gain = std::make_unique<MfemMeshFunction>(space, gainFile);
            xRayGain = std::make_unique<XRayGain>(*gain);
            exchangeController.addModel(xRayGain.get());
        }

        std::cout << stringifyAbsorptionSummary(
                exchangeController.absorb(intersections, generateInitialEnergies(laser), absorbedEnergy)
        );
    }

    std::ofstream raysFile(raysOutputFilename);
    raysFile << stringifyRaysToMsgpack(allIntersections);
    std::ofstream absorbedResult(absorbedEnergyFilename);
    absorbedResult << absorbedEnergy;
}