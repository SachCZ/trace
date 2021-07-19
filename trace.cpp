#include <raytracer.h>
#include "yaml-cpp/yaml.h"
#include "trace_config.h"
#include <fstream>

void divideByVolume(const raytracer::MfemMesh& mesh, raytracer::MfemMeshFunction& func){
    for (auto element : mesh.getElements()){
        auto volume = raytracer::getElementVolume(*element);
        func.setValue(*element, func.getValue(*element) / volume);
    }
}

int main(int argc, char *argv[]) {
    using namespace raytracer;

    for (int argNumber = 1; argNumber < argc; argNumber++){
        YAML::Node config;
        try {
            config = YAML::LoadFile(argv[argNumber]);
        } catch (const YAML::BadFile& err){
            std::cerr << "File " <<  argv[argNumber] << " not found." << std::endl;
            return 1;
        }
        TraceConfig traceConfig(config);

        IntersectionSet allIntersections;
        MfemMeshFunction absorbedPower(*traceConfig.l2Space);

        for (const auto &laserConfig : traceConfig.laserConfigs) {
            auto initialDirections = generateInitialDirections(laserConfig.laser);
            auto intersections = findIntersections(
                    *traceConfig.mesh,
                    initialDirections,
                    {*laserConfig.totalReflect, *laserConfig.snellsLaw},
                    intersectStraight,
                    dontStop
            );
            allIntersections.insert(allIntersections.end(), intersections.begin(), intersections.end());

            auto initialPowers = generateInitialPowers(laserConfig.laser);
            auto modelPowers = laserConfig.exchangeController.genPowers(intersections, initialPowers);
            auto rayPowers = modelPowersToRayPowers(modelPowers, initialPowers);
            absorbRayPowers(absorbedPower, rayPowers, intersections);
            if (!laserConfig.energiesOutputFilename.empty()){
                std::ofstream energiesFile(laserConfig.energiesOutputFilename);
                rayPowersToMsgpack(rayPowers, energiesFile);
            }
        }

        if (!traceConfig.raysOutputFilename.empty()){
            std::ofstream raysFile(traceConfig.raysOutputFilename);
            raysFile << stringifyRaysToMsgpack(allIntersections);
        }

        if (!traceConfig.absorbedPowerFilename.empty()){
            std::ofstream absorbedResult(traceConfig.absorbedPowerFilename);
            absorbedResult << absorbedPower;
        }

        if (!traceConfig.meshFilename.empty()){
            std::ofstream meshResult(traceConfig.meshFilename);
            traceConfig.mesh->getMfemMesh()->Print(meshResult);
        }
    }
}