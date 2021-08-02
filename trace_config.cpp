#include "trace_config.h"
#include "expression.cpp" //hack to lower compile time, not a big fan
#include <random>

LaserTraceConfig::LaserTraceConfig(const YAML::Node &laserConfig, const TraceConfig &traceConfig) {
    using namespace raytracer;
    gradient = parseGradient(laserConfig, traceConfig);
    laser = parseLaser(laserConfig);
    reflectedMarker = std::make_unique<Marker>();
    criticalMarker = std::make_unique<Marker>();

    if (traceConfig.temperature && traceConfig.ionization) {
        frequency = std::make_unique<MfemMeshFunction>(*traceConfig.l2Space, [&](const Element &e) {
            return calcSpitzerFreq(
                    traceConfig.electronDensity->getValue(e),
                    traceConfig.temperature->getValue(e),
                    traceConfig.ionization->getValue(e),
                    this->laser.wavelength
            );
        });
    } else if (!traceConfig.temperature && !traceConfig.ionization) {
        frequency = std::make_unique<MfemMeshFunction>(*traceConfig.l2Space, [&](const Element &) { return 0; });
    }
    refractiveIndex = std::make_unique<MfemMeshFunction>(*traceConfig.l2Space, [&](const Element &e) {
        return calcRefractIndex(
                traceConfig.electronDensity->getValue(e),
                laser.wavelength,
                frequency->getValue(e)
        );
    });

    snellsLaw = std::make_unique<SnellsLawBend>(traceConfig.mesh.get(), refractiveIndex.get(), gradient.get());
    totalReflect = std::make_unique<TotalReflect>(traceConfig.mesh.get(), refractiveIndex.get(), gradient.get(),
                                                  reflectedMarker.get());
    reflectOnCritical = std::make_unique<ReflectOnCritical>(
            traceConfig.mesh.get(),
            refractiveIndex.get(),
            traceConfig.electronDensity.get(),
            calcCritDens(laser.wavelength).asDouble,
            gradient.get(),
            criticalMarker.get()
    );
    exchangeController = parseExchangeController(laserConfig, traceConfig);
    if (laserConfig["energies_filename"])
        energiesOutputFilename = parse<std::string>(laserConfig, "energies_filename");
}

raytracer::PowerExchangeController
LaserTraceConfig::parseExchangeController(const YAML::Node &laserConfig, const TraceConfig &traceConfig) {
    using namespace raytracer;
    PowerExchangeController result;
    if (!laserConfig["power_exchange"]) {
        zeroExchange = std::make_unique<ZeroExchange>();
        result.addModel(zeroExchange.get());
        return result;
    }
    for (const auto &node: laserConfig["power_exchange"]) {
        auto type = node.as<std::string>();
        if (type == "gain") {
            std::ifstream gainFile(parse<std::string>(laserConfig, "gain_file"));
            gain = std::make_unique<MfemMeshFunction>(*traceConfig.l2Space, gainFile);
            xRayGain = std::make_unique<XRayGain>(*gain);
            result.addModel(xRayGain.get());
        } else if (type == "bremsstrahlung") {
            invBremssCoeff = std::make_unique<MfemMeshFunction>(*traceConfig.l2Space, [&](const Element &e) {
                return calcInvBremssCoeff(
                        traceConfig.electronDensity->getValue(e),
                        laser.wavelength,
                        frequency->getValue(e)
                );
            });

            bremsstrahlung = std::make_unique<Bremsstrahlung>(invBremssCoeff.get());
            result.addModel(bremsstrahlung.get());
        } else if (type == "resonance") {
            resonance = std::make_unique<Resonance>(laser.wavelength, reflectedMarker.get(), gradient.get());
            result.addModel(resonance.get());
        }
    }
    return result;
}

raytracer::Laser LaserTraceConfig::parseLaser(const YAML::Node &laserConfig) {
    using namespace raytracer;
    Length wavelength{parse<double>(laserConfig, "wavelength")};
    auto direction{parseXY<Vector>(laserConfig, "direction")};
    auto startPoint = parseXY<Point>(laserConfig, "start_point");
    auto endPoint = parseXY<Point>(laserConfig, "end_point");
    auto raysCount = parse<int>(laserConfig, "rays_count");
    Laser::PowerFun powerFunction;
    if (laserConfig["power"] && laserConfig["spatial_FWHM"]) {
        auto spatialFWHM = parse<double>(laserConfig, "spatial_FWHM");
        Power power{parse<double>(laserConfig, "power")};
        powerFunction = raytracer::MaxValGaussian(spatialFWHM, power.asDouble, 0);
    } else if (laserConfig["constant_power"]) {
        auto value = parse<double>(laserConfig, "constant_power");
        powerFunction = [value](double) {return value;};
    } else {
        powerFunction = [](double) { return 0; };
    }
    return raytracer::Laser{
            wavelength,
            [direction](const Point &) { return direction; },
            powerFunction,
            startPoint,
            endPoint,
            raysCount
    };
}

void vec_bdr(const mfem::Vector &point, mfem::Vector &result) {
    result[0] = 1.7145e+24 * point[0];
    result[1] = 0;
}

std::unique_ptr<raytracer::Gradient>
LaserTraceConfig::parseGradient(const YAML::Node &node, const TraceConfig &traceConfig) {
    using namespace raytracer;
    auto gradType = parse<std::string>(node, "grad_type");

    raytracer::VectorField gradAtPoints;
    if (gradType == "ls") {
        gradAtPoints = calcHousGrad(*traceConfig.mesh, *traceConfig.electronDensity, false);
    } else if (gradType == "mfem") {
        auto bdrXString = parse<std::string>(node, "grad_boundary_x");
        auto bdrYString = parse<std::string>(node, "grad_boundary_y");
        Expression expressionX(bdrXString);
        Expression expressionY(bdrYString);
        mfem::VectorFunctionCoefficient gradientBoundaryValue(2, vec_bdr);
        /**
        [&](const mfem::Vector &point,
               mfem::Vector &result) {
result[0] = expressionX(point[0], point[1]);
result[1] = expressionY(point[0], point[1]);
});
         */
        gradAtPoints = mfemGradient(*traceConfig.mesh, *traceConfig.electronDensity, &gradientBoundaryValue);
    } else if (gradType == "integral") {
        gradAtPoints = calcIntegralGrad(*traceConfig.mesh, *traceConfig.electronDensity);
    } else {
        throw ParsingError("Unknown gradient type");
    }
    gradAtPoints = raytracer::setValue(
            gradAtPoints,
            traceConfig.mesh->getBoundaryPoints(),
            {1, 0}
    );
    return std::make_unique<LinInterGrad>(gradAtPoints);
}

TraceConfig::TraceConfig(const YAML::Node &config) {
    mesh = parseMesh(config);
    l2Space = std::make_unique<raytracer::MfemL20Space>(*this->mesh);
    ionization = parseIonization(config);
    electronDensity = parseDensity(config);
    temperature = parseTemperature(config);
    laserConfigs = parseLasers(config);
    if (config["rays_filename"])
        raysOutputFilename = parse<std::string>(config, "rays_filename");
    if (config["absorbed_power_filename"])
        absorbedPowerFilename = parse<std::string>(config, "absorbed_power_filename");
    if (config["mesh_output_filename"])
        meshFilename = parse<std::string>(config, "mesh_output_filename");
}

std::vector<LaserTraceConfig> TraceConfig::parseLasers(const YAML::Node &config) {
    if (!config["lasers"]) throw ParsingError("Lasers not specified");
    std::vector<LaserTraceConfig> result;
    for (const auto &node : config["lasers"]) {
        result.emplace_back(LaserTraceConfig{node, *this});
    }
    return result;
}

void randomizeMesh(raytracer::MfemMesh &mesh, double factor, int xSegments, int ySegments) {
    raytracer::MfemMesh::Displacements displacements;
    std::mt19937 gen(std::random_device{}());
    double maxXDisplacement = 1.0 / xSegments * factor;
    double maxYDisplacement = 1.0 / ySegments * factor;
    std::uniform_real_distribution xDist(-maxXDisplacement, maxXDisplacement);
    std::uniform_real_distribution yDist(-maxYDisplacement, maxYDisplacement);

    const auto &innerPoints = mesh.getInnerPoints();
    for (auto point : mesh.getPoints()) {
        if (std::find(innerPoints.begin(), innerPoints.end(), point) == innerPoints.end()) {
            displacements.emplace_back(raytracer::Vector{0, 0});
        } else {
            displacements.emplace_back(raytracer::Vector{xDist(gen), yDist(gen)});
        }
    }
    mesh.moveNodes(displacements);
}

TraceConfig::MeshPtr TraceConfig::parseMesh(const YAML::Node &config) {
    if (config["mesh_file"]) {
        auto meshFilename = config["mesh_file"].as<std::string>();
        return std::make_unique<raytracer::MfemMesh>(meshFilename);
    } else if (config["mesh"]) {
        if (
                !config["mesh"]["x0"] ||
                !config["mesh"]["x1"] ||
                !config["mesh"]["y0"] ||
                !config["mesh"]["y1"] ||
                !config["mesh"]["x_segments"] ||
                !config["mesh"]["y_segments"]

                )
            throw ParsingError("Invalid mesh config");

        auto x0 = config["mesh"]["x0"].as<double>();
        auto x1 = config["mesh"]["x1"].as<double>();
        auto y0 = config["mesh"]["y0"].as<double>();
        auto y1 = config["mesh"]["y1"].as<double>();
        auto xSegments = config["mesh"]["x_segments"].as<size_t>();
        auto ySegments = config["mesh"]["y_segments"].as<size_t>();

        raytracer::SegmentedLine sideX{x0, x1, xSegments};
        raytracer::SegmentedLine sideY{y0, y1, ySegments};
        auto result = std::make_unique<raytracer::MfemMesh>(sideX, sideY);

        if (config["mesh"]["randomize_factor"]) {
            auto randomFactor = config["mesh"]["randomize_factor"].as<double>();
            randomizeMesh(*result, randomFactor, xSegments, ySegments);
            result->updateMesh();
        }
        return result;
    } else {
        throw ParsingError("Mesh not specified");
    }
}

FunctionPtr TraceConfig::parseIonization(const YAML::Node &config) const {
    if (config["ionization_file"]) {
        return parseFunction(config, "ionization_file");
    } else if (config["ionization_profile"]) {
        auto expression = Expression(config["ionization_profile"].as<std::string>());
        return std::make_unique<raytracer::MfemMeshFunction>(*l2Space, [&](const raytracer::Point &point) {
            return expression(point.x, point.y);
        });
    } else {
        return nullptr;
    }
}

FunctionPtr TraceConfig::parseTemperature(const YAML::Node &config) const {
    if (config["temperature_file"]) {
        return parseFunction(config, "temperature_file");
    } else if (config["temperature_profile"]) {
        auto expression = Expression(config["temperature_profile"].as<std::string>());
        return std::make_unique<raytracer::MfemMeshFunction>(*l2Space, [&](const raytracer::Point &point) {
            return expression(point.x, point.y);
        });
    } else {
        return nullptr;
    }
}


FunctionPtr TraceConfig::parseDensity(const YAML::Node &config) const {
    if (config["ele_dens_file"]) {
        return parseFunction(config, "ele_dens_file");
    } else if (config["ele_dens_profile"]) {
        auto expression = Expression(config["ele_dens_profile"].as<std::string>());
        return std::make_unique<raytracer::MfemMeshFunction>(*l2Space, [&](const raytracer::Point &point) {
            return expression(point.x, point.y);
        });
    } else if ((config["ion_dens_file"] || config["ion_dens_profile"]) && config["atomic_mass"]) {
        FunctionPtr density;
        if (config["ion_dens_file"]) {
            density = parseFunction(config, "ion_dens_file");
        } else if (config["ion_dens_profile"]) {
            auto expression = Expression(config["ion_dens_profile"].as<std::string>());
            density = std::make_unique<raytracer::MfemMeshFunction>(*l2Space, [&](const raytracer::Point &point) {
                return expression(point.x, point.y);
            });
        }
        auto A = config["atomic_mass"].as<double>();
        return std::make_unique<raytracer::MfemMeshFunction>(*l2Space, [&](const raytracer::Element &e) {
            double massUnit = 1.6605e-24;
            return density->getValue(e) * ionization->getValue(e) / massUnit / A;
        });
    } else {
        throw ParsingError("Density not specified use: ele_dens_file, ele_dens_profile,"
                           " ion_dens_file with atomic_mass or ion_dens_profile with atomic_mass.");
    }
}

FunctionPtr TraceConfig::parseFunction(const YAML::Node &config, const std::string &configFunctionName) const {
    auto filename = config[configFunctionName].as<std::string>();
    std::ifstream file(filename);
    return std::make_unique<raytracer::MfemMeshFunction>(*l2Space, file);
}
