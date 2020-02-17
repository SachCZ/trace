#include <raytracer/geometry/Mesh.h>
#include <raytracer/physics/MathFunctions.h>
#include <raytracer/physics/Laser.h>
#include <raytracer/geometry/GeometryFunctions.h>
#include <raytracer/geometry/Ray.h>
#include <raytracer/geometry/MeshFunction.h>
#include <chrono>
#include "mfem.hpp"

using namespace raytracer::geometry;
using namespace raytracer::physics;

struct StopAtCritical {
    explicit StopAtCritical(const MeshFunction &density) :
            density(density) {}

    bool operator()(const Intersection &intersection, const LaserRay &laserRay) {
        const auto element = *intersection.nextElement;

        auto currentDensity = density.getValue(element);
        auto criticalDensity = laserRay.getCriticalDensity();
        return currentDensity > criticalDensity.asDouble;
    }

private:
    const MeshFunction &density;
};


std::unique_ptr<Intersection> continueStraight(const Intersection &intersection) {
    return findClosestIntersection(intersection.orientation, intersection.nextElement->getFaces(), intersection.face);
}

class GradientCalculator {
public:
    virtual Vector getGradient(const MeshFunction &, const Point &) const = 0;
};

class ConstantGradientCalculator : public GradientCalculator {
public:
    explicit ConstantGradientCalculator(const Vector &gradient) : gradient(gradient) {}

    Vector getGradient(const MeshFunction &, const Point &) const override {
        return this->gradient;
    }

private:
    const Vector gradient;
};

struct SnellsLaw {
    explicit SnellsLaw(const MeshFunction &density, const GradientCalculator &gradientCalculator) :
            density(density), gradientCalculator(gradientCalculator) {}

    std::unique_ptr<Intersection> operator()(const Intersection &intersection, const LaserRay &laserRay) {
        const auto previousElement = intersection.previousElement;
        const auto nextElement = intersection.nextElement;

        if (!previousElement) {
            return findClosestIntersection(intersection.orientation, nextElement->getFaces(), intersection.face);
        }

        auto newIntersection = findClosestIntersection(
                intersection.orientation,
                nextElement->getFaces(),
                intersection.face);
        newIntersection->orientation.direction = getDirection(intersection, laserRay.getCriticalDensity());
        return newIntersection;
    }

private:
    const MeshFunction &density;
    const GradientCalculator &gradientCalculator;

    Vector getDirection(const Intersection &intersection, Density criticalDensity) {
        const double n1 = std::sqrt(1 - density.getValue(*intersection.previousElement) / criticalDensity.asDouble);
        const double n2 = std::sqrt(1 - density.getValue(*intersection.nextElement) / criticalDensity.asDouble);

        const auto gradient = gradientCalculator.getGradient(density, intersection.orientation.point);
        const auto &direction = intersection.orientation.direction;

        const auto l = 1 / direction.getNorm() * direction;
        auto n = 1 / gradient.getNorm() * gradient;
        auto c = (-1) * n * l;
        if (c < 0) {
            c = -c;
            n = (-1) * n;
        }
        const double r = n1 / n2;

        auto root = 1 - r * r * (1 - c * c);
        if (root > 0) {
            return r * l + (r * c - std::sqrt(root)) * n;
        } else {
            return l + 2 * c * n;
        }
    }
};


double density(const mfem::Vector &x) {
    return 12.8e20 * x(0);
}

struct EndAbsorber {
    explicit EndAbsorber(MeshFunction &absorbedEnergy) : absorbedEnergy(absorbedEnergy) {}

    void operator()(const LaserRay &laserRay) {
        const auto &intersection = laserRay.intersections.back();
        const auto element = intersection.previousElement;
        if (!element) return;
        absorbedEnergy.addValue(*element, laserRay.energy.asDouble);
    }

private:
    MeshFunction &absorbedEnergy;
};

int main(int, char *[]) {


    //MFEM boilerplate -------------------------------------------------------------------------------------------------
    //DiscreteLine side{};
    //side.segmentCount = 100;
    //side.length = 1;
    //auto mfemMesh = constructRectangleMesh(side, side);

    auto mfemMesh = std::make_unique<mfem::Mesh>("test_mesh.vtk", 1, 0);


    mfem::L2_FECollection finiteElementCollection(0, 2);
    mfem::FiniteElementSpace finiteElementSpace(mfemMesh.get(), &finiteElementCollection);

    mfem::GridFunction densityGridFunction(&finiteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    mfem::GridFunction absorbedEnergyGridFunction(&finiteElementSpace);
    absorbedEnergyGridFunction = 0;
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());
    MeshFunction densityMeshFunction(densityGridFunction, finiteElementSpace);
    MeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, finiteElementSpace);

    ConstantGradientCalculator gradientCalculator(Vector(12.8e20, 0));

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point) { return Vector(1, -0.3); },
            Gaussian(0.1),
            Point(-1.1, 0.9),
            Point(-1.1, 0.7)
    );

    laser.generateRays(1);
    laser.generateIntersections(
            mesh, SnellsLaw(densityMeshFunction, gradientCalculator),
            StopAtCritical(densityMeshFunction));

    EndAbsorber endAbsorber(absorbedEnergyMeshFunction);
    for (const auto &laserRay : laser.getRays()) {
        endAbsorber(laserRay);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto nanosecondDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    auto duration = nanosecondDuration * 1e-9;

    std::cout << "Execution took: " << duration << "s." << std::endl;

    laser.saveRaysToJson("rays.json");

    //GLVIS-------------------------------------------------------------------------------------------------------------
    char vishost[] = "localhost";
    int visport = 19916;
    mfem::socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
}