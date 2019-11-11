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
    Mesh mesh("./mesh.stl");
    Ray ray(Point(-0.995, 5), Vector(0, -1));
    ray.traceThrough(mesh, [](const Intersection& intersection){
        //auto x = intersection.point.x;
        //auto y = intersection.point.y;
        return Vector(0, -1);
    });
    ray.saveToJson("ray");
    mesh.saveToJson("mesh");
}