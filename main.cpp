#include "raytracer/geometry/Mesh.h"
#include "raytracer/geometry/Ray.h"
#include "raytracer/geometry/Point.h"
#include "raytracer/geometry/Vector.h"
#include <cmath>

using Mesh = raytracer::geometry::Mesh;
using Ray = raytracer::geometry::Ray;
using Point = raytracer::geometry::Point;
using Vector = raytracer::geometry::Vector;
using RayState = raytracer::geometry::RayState;

int main(int argc, char* argv[])
{
    Mesh mesh("./broken_mesh.vtk");
    Ray ray(Point(-2, 0.9), Vector(1, 0));
    ray.traceThrough(mesh, [](const RayState & rayState){
        auto point = rayState.lastIntersection.point;
        auto x = point.x;
        auto y = point.y;
        return Vector(3*y, -1);
    }, [](const RayState& rayState){return false;});
    ray.saveToJson("ray");
    mesh.saveToJson("mesh");
}