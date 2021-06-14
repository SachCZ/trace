#include <raytracer.h>
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <random>
#include "exprtk.hpp"

using namespace raytracer;

VectorField
calcAnalyticGrad(const VectorField &discreteGradient, const std::function<Vector(Point)> &analyticGradFunc) {
    VectorField result;
    for (const auto &pointAndVector : discreteGradient) {
        auto point = pointAndVector.first;
        result[point] = analyticGradFunc(*point);
    }
    return result;
}

std::unique_ptr<std::ostream> openFile(const std::vector<std::string> &paths) {
    std::stringstream outputFilename;
    for (const auto &path : paths) {
        outputFilename << path;
    }
    return std::make_unique<std::ofstream>(outputFilename.str());
}

void writeGradResult(
        const std::string &basePath,
        const std::string &fileId,
        const MfemMesh &mesh,
        const MfemMeshFunction &func,
        const VectorField &gradient,
        const VectorField &analyticGradient
) {
    using namespace raytracer;

    *openFile({basePath, "mesh", fileId, ".mfem"}) << mesh;
    *openFile({basePath, "func", fileId, ".gf"}) << func;
    *openFile({basePath, "grad", fileId, ".msgpack"}) << gradient;
    *openFile({basePath, "analytic_grad", fileId, ".msgpack"}) << analyticGradient;

    auto dualMeshOutput = openFile({basePath, "dual_mesh", fileId, ".mfem"});
    writeDualMesh(*dualMeshOutput, mesh);
}

int main(int, char *argv[]) {
    YAML::Node config = YAML::LoadFile(argv[1]);

    std::function<double(Point)> analyticFunc;
    std::function<Vector(Point)> analyticGradFunc;
    if (config["function"]) {
        auto functionName = config["function"].as<std::string>();
        if (functionName == "sin") {
            using namespace std;
            analyticFunc = [](const Point &point) {
                return 1.62 + 0.31 * sin(M_PI * (3.79 * point.x + 2.98 * point.y));
            };
            analyticGradFunc = [](const Point &point) {
                auto x = point.x;
                auto y = point.y;
                return Vector{3.69106 * cos(11.9066 * x + 9.36195 * y), 2.9022 * cos(11.9066 * x + 9.36195 * y)};
            };
        } else if (functionName == "lin") {
            using namespace std;
            analyticFunc = [](const Point &point) {
                return 2 * point.x + 3 * point.y;
            };
            analyticGradFunc = [](const Point &) {
                return Vector{2, 3};
            };
        }
    } else {
        throw std::logic_error("'function' not defined in config file");
    }

    if (
            config["meshes"]["segments"] &&
            config["meshes"]["random_factors"]
            ) {
        auto type = config["type"].as<std::string>();
        for (const auto &segmentsNode : config["meshes"]["segments"]) {
            auto count = segmentsNode.as<int>();
            std::stringstream basePath;
            auto mfemMesh = std::make_unique<mfem::Mesh>(
                    count, count,
                    mfem::Element::QUADRILATERAL,
                    true,
                    1.0, 1.0,
                    true
            );
            auto verticesCount = mfemMesh->GetNV();
            mfem::Vector displacements(verticesCount * 2);
            std::mt19937 gen(std::random_device{}());
            auto randomFactors = config["meshes"]["random_factors"];
            for (auto randomFactorNode : randomFactors) {
                auto randomFactor = randomFactorNode.as<double>();
                double maxDisplacement = 1.0 / count * randomFactor;
                std::uniform_real_distribution dist(-maxDisplacement, maxDisplacement);
                for (int i = 0; i < verticesCount * 2; i++) {
                    if (
                            i % (count + 1) == 0 ||
                            (i + 1) % (count + 1) == 0 ||
                            i < count + 1 ||
                            (i > verticesCount - (count + 1) && i < verticesCount + (count + 1)) ||
                            i > 2 * verticesCount - (count + 1)
                            ) {
                        displacements[i] = 0;
                    } else {
                        displacements[i] = dist(gen);
                    }
                }
                mfemMesh->MoveVertices(displacements);
                MfemMesh mesh(mfemMesh.get());

                MfemL20Space space(mesh);
                MfemMeshFunction func(space, analyticFunc);
                VectorField gradient;
                if (type == "haus") {
                    gradient = calcHousGrad(mesh, func);
                }
                if (type == "integral") {
                    gradient = calcIntegralGrad(mesh, func);
                }
                if (type == "mfem") {
                    mfem::H1_FECollection h1FiniteElementCollection{1, 2};
                    mfem::FiniteElementSpace h1FiniteElementSpace(mesh.getMfemMesh(), &h1FiniteElementCollection, 2);
                    mfem::FunctionCoefficient boundaryValue([&](const mfem::Vector &point) {
                        return analyticFunc({point[0], point[1]});
                    });
                    mfem::VectorFunctionCoefficient gradientBoundaryValue(2, [&](const mfem::Vector &point,
                                                                                 mfem::Vector &result) {
                        auto value = analyticGradFunc({point[0], point[1]});
                        result[0] = value.x;
                        result[1] = value.y;
                    });
                    gradient = mfemGradient(
                            mesh,
                            func,
                            &gradientBoundaryValue,
                    0,
                            1.0 / count
                    );
                }
                auto analyticGradient = calcAnalyticGrad(gradient, analyticGradFunc);
                writeGradResult("output/", std::to_string(count) + "_" + randomFactorNode.as<std::string>(), mesh, func,
                                gradient, analyticGradient);
            }
        }
    } else {
        throw std::logic_error("meshes have invalid config in config file");
    }
}


