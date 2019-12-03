#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/Functions.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/Functions.h>
#include <raytracer/geometry/Ray.h>
#include "mfem.hpp"

using namespace raytracer::geometry;
using namespace raytracer::physics;

double density(const mfem::Vector &x){
    return 12.8e20 * x(0);
}


int main(int argc, char *argv[]) {

    DiscreteLine side{};
    side.segmentCount = 100;
    side.width = 1;
    Mesh mesh(side, side);

    mfem::Mesh externalMesh(
            side.segmentCount,
            side.segmentCount,
            mfem::Element::Type::QUADRILATERAL,
            true,
            side.width,
            side.width,
            true);

    mfem::H1_FECollection finiteElementCollection(1, 2);
    mfem::FiniteElementSpace finiteElementSpace(&externalMesh, &finiteElementCollection);
    mfem::GridFunction densityGridFunction(&finiteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    Laser laser(Length{1315e-7}, [](const Point &point) { return Vector(1, 0); }, Gaussian(0.2));

    std::vector<std::vector<Intersection>> rays;
    for (const auto &laserRay : laser.generateRays(100, Point(-1, 0.1), Point(-1, 0.9))) {
        Ray ray(HalfLine{laserRay.startPoint, laserRay.direction});

        auto intersections = ray.findIntersections(
                mesh,
                [](const Intersection& intersection, const Element& element) {
                    return findClosestIntersection(intersection.orientation, element.getFaces());
                },
                [&laserRay, &densityGridFunction](const Intersection& intersection, const Element& element) {
                    mfem::Array<double> nodalValues;
                    mfem::IntegrationPoint integrationPoint{};
                    integrationPoint.Set2(intersection.orientation.point.x, intersection.orientation.point.y);
                    auto currentDensity = densityGridFunction.GetValue(element.getId(), integrationPoint);
                    auto criticalDensity = laserRay.getCriticalDensity().asDouble;
                    return currentDensity > criticalDensity;
                }
        );
        rays.emplace_back(intersections);
    }
    std::cout << "";
}