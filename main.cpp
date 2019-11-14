#include "raytracer/geometry/Mesh.h"
#include "raytracer/geometry/Ray.h"
#include "raytracer/geometry/Point.h"
#include "raytracer/geometry/Vector.h"
#include <cmath>

using Mesh = raytracer::geometry::Mesh;
using Ray = raytracer::geometry::Ray;
using Point = raytracer::geometry::Point;
using Vector = raytracer::geometry::Vector;
using Intersection = raytracer::geometry::Intersection;

int main(int argc, char* argv[])
{
    Mesh mesh("./broken_mesh.vtk");
    Ray ray(Point(-2, 0.9), Vector(1, 0));
    ray.traceThrough(mesh, [](const Intersection& intersection){
        auto x = intersection.point.x;
        auto y = intersection.point.y;
        return Vector(3*y, -1);
    });
    ray.saveToJson("ray");
    mesh.saveToJson("mesh");
}