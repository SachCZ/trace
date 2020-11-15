#include <raytracer.h>

void gradientOnMesh(size_t segmentsCount, const std::function<double(raytracer::Point)>& analytic){
    using namespace raytracer;

    MfemMesh mesh(SegmentedLine{1.0, segmentsCount}, SegmentedLine{1.0, segmentsCount});
    MfemL20Space space(mesh);
    MfemMeshFunction func(space, analytic);
    std::stringstream outputFilename;
    outputFilename << "output/household" << segmentsCount << ".msgpack";
    std::ofstream output(outputFilename.str());
    auto gradient = calcHousGrad(mesh, func);
    output << gradient;
    std::stringstream functionOutputFilename;
    functionOutputFilename << "output/household" << segmentsCount << ".gf";
    std::ofstream functionOutput(functionOutputFilename.str());
    functionOutput << func;

    std::stringstream meshOutputFilename;
    meshOutputFilename << "output/household" << segmentsCount << ".mfem";
    std::ofstream meshOutput(meshOutputFilename.str());
    meshOutput << mesh;
}

int main(int, char *[]) {
    for (size_t count = 10; count <= 100; count += 10){
        gradientOnMesh(count, [](const raytracer::Point& point){
            return 1.62 + 0.31*std::sin(M_PI*(3.79*point.x + 2.98*point.y));
        });
    }
}