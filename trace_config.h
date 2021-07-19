#ifndef TRACE_TRACE_CONFIG_H
#define TRACE_TRACE_CONFIG_H

#include <yaml-cpp/yaml.h>
#include <raytracer.h>


using FunctionPtr = std::unique_ptr<raytracer::MfemMeshFunction>;

template<typename T>
static T parse(const YAML::Node &node, const std::string &property) {
    if (!node[property]) throw std::logic_error(property + " is required");
    return node[property].as<T>();
}

struct ParsingError : public std::exception
{
    explicit ParsingError(const char * msg): msg(msg) {}

    [[nodiscard]] const char * what () const noexcept override
    {
        return msg;
    }

private:
    const char * msg;
};

class TraceConfig;

class LaserTraceConfig {
public:
    LaserTraceConfig(const YAML::Node &laserConfig, const TraceConfig &traceConfig);

    std::unique_ptr<raytracer::Gradient> gradient;
    raytracer::Laser laser;
    FunctionPtr frequency;
    FunctionPtr refractiveIndex;
    raytracer::PowerExchangeController exchangeController;
    std::unique_ptr<raytracer::SnellsLawBend> snellsLaw;
    std::unique_ptr<raytracer::TotalReflect> totalReflect;
    std::unique_ptr<raytracer::ReflectOnCritical> reflectOnCritical;
    std::string energiesOutputFilename{};

private:
    FunctionPtr gain;
    FunctionPtr invBremssCoeff;
    std::unique_ptr<raytracer::Marker> reflectedMarker;
    std::unique_ptr<raytracer::Marker> criticalMarker;

    //Exchange models
    std::unique_ptr<raytracer::XRayGain> xRayGain;
    std::unique_ptr<raytracer::Bremsstrahlung> bremsstrahlung;
    std::unique_ptr<raytracer::Resonance> resonance;
    std::unique_ptr<raytracer::ZeroExchange> zeroExchange;

    raytracer::PowerExchangeController parseExchangeController(
            const YAML::Node &laserConfig,
            const TraceConfig &traceConfig
    );

    raytracer::Laser parseLaser(const YAML::Node &laserConfig);

    template<typename T>
    T parseXY(const YAML::Node &node, const std::string &property) {
        if (!node[property]["x"] || !node[property]["y"]) {
            throw std::logic_error(property + " not specified or ill formatted, use x: y:");
        }
        return {node[property]["x"].as<double>(), node[property]["y"].as<double>()};
    }

    static std::unique_ptr<raytracer::Gradient> parseGradient(const YAML::Node &node, const TraceConfig &traceConfig);

};

class TraceConfig {
    using MeshPtr = std::unique_ptr<raytracer::MfemMesh>;
    using SpacePtr = std::unique_ptr<raytracer::MfemL20Space>;
public:
    explicit TraceConfig(const YAML::Node &config);

    MeshPtr mesh;
    SpacePtr l2Space;
    FunctionPtr ionization;
    FunctionPtr electronDensity;
    FunctionPtr temperature;
    std::vector<LaserTraceConfig> laserConfigs;
    std::string raysOutputFilename{};
    std::string absorbedPowerFilename{};
    std::string meshFilename{};

private:
    std::vector<LaserTraceConfig> parseLasers(const YAML::Node &config);

    static MeshPtr parseMesh(const YAML::Node &config);

    [[nodiscard]] FunctionPtr parseIonization(const YAML::Node &config) const;

    [[nodiscard]] FunctionPtr parseTemperature(const YAML::Node &config) const;

    [[nodiscard]] FunctionPtr parseDensity(const YAML::Node &config) const;

    [[nodiscard]] FunctionPtr parseFunction(const YAML::Node &config, const std::string &configFunctionName) const;
};

#endif //TRACE_TRACE_CONFIG_H
