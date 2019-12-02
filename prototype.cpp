#include <utility>

#include "raytracer/geometry/Mesh.h"
#include "raytracer/geometry/MeshFunction.h"
#include "raytracer/geometry/Ray.h"
#include "raytracer/geometry/Point.h"
#include "raytracer/geometry/Vector.h"
#include "raytracer/physics/LaserRay.h"
#include "raytracer/physics/Magnitudes.h"
#include "raytracer/physics/Functions.h"
#include "raytracer/physics/Laser.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <functional>

using Mesh = raytracer::geometry::Mesh;
using Ray = raytracer::geometry::Ray;
using Point = raytracer::geometry::Point;
using Vector = raytracer::geometry::Vector;
using Quadrilateral = raytracer::geometry::Quadrilateral;
using RayState = raytracer::geometry::RayState;
using LaserRay = raytracer::physics::LaserRay;
using Laser = raytracer::physics::Laser;
using MeshFunction = raytracer::geometry::MeshFunction;
using Energy = raytracer::physics::Energy;
using Length = raytracer::physics::Length;
using LinearFunction = raytracer::physics::LinearFunction;
using Gaussian = raytracer::physics::Gaussian;
namespace JSON = raytracer::utility::json;
namespace utility = raytracer::utility;

int main(int argc, char *argv[]) {
    Mesh mesh("./broken_mesh.vtk");

    MeshFunction absorbedEnergy(mesh);
    MeshFunction density(mesh);

    LinearFunction linearFunction(6.3e20, 6.4e20);
    density.setAll(mesh.getQuads(), [&](const Quadrilateral &quad) { return linearFunction(quad.getAveragePoint().x);});

    Laser laser(Length{1315e-7}, [](const Point& point){ return Vector(1, 0);}, Gaussian(0.2));

    std::ofstream file("ray.json");
    JSON::Value rays(JSON::arrayValue);

    for (const auto &laserRay : laser.generateRays(100, Point(-2, -0.9), Point(-2, 0.9))) {
        Ray ray(laserRay.startPoint, laserRay.direction);

        auto nextDirection = [](const RayState &rayState) {
            return rayState.currentDirection;
        };

        auto stopCondition = density.greaterOrEqual(
                laserRay.getCriticalDensity().asDouble,
                [&](const Quadrilateral& quad) { absorbedEnergy[quad] += laserRay.energy.asDouble; });

        ray.traceThrough(mesh, nextDirection, stopCondition);

        rays.append(ray.getJsonValue());
    }
    JSON::Value root;
    root["rays"] = rays;
    root["absorbedEnergy"] = absorbedEnergy.getJsonValue();
    file << root;
    mesh.saveToJson("mesh");
}