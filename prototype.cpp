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

struct DontStop {
    bool operator()(const Intersection&, const LaserRay&) {
        return false;
    }
};


std::unique_ptr<Intersection> continueStraight(const Intersection &intersection, const LaserRay &) {
    return findClosestIntersection(intersection.orientation, intersection.nextElement->getFaces(), intersection.face);
}

class GradientCalculator {
public:
    virtual Vector getGradient(const Intersection &) const = 0;
};

class ConstantGradientCalculator : public GradientCalculator {
public:
    explicit ConstantGradientCalculator(const Vector &gradient) : gradient(gradient) {}

    Vector getGradient(const Intersection &) const override {
        return this->gradient;
    }

private:
    const Vector gradient;
};

class H1GradientCalculator : public GradientCalculator {
public:
    H1GradientCalculator(mfem::FiniteElementSpace &l2Space, mfem::FiniteElementSpace &h1Space):
            l2Space(l2Space), h1Space(h1Space), _density(&h1Space) {}

    Vector getGradient(const Intersection &intersection) const override {
        auto point = intersection.orientation.point;
        auto previousGradient = this->getGradientAt(*intersection.previousElement, point);
        auto nextGradient = this->getGradientAt(*intersection.nextElement, point);
        return 0.5 * (previousGradient + nextGradient);
    }

    void updateDensity(mfem::GridFunction& density){
        this->_density = convertH1toL2(density);
    }

private:
    mfem::FiniteElementSpace &l2Space;
    mfem::FiniteElementSpace &h1Space;
    mfem::GridFunction _density;

    Vector getGradientAt(const Element& element, const Point& point) const {
        mfem::Vector result(2);
        mfem::IntegrationPoint integrationPoint{};
        integrationPoint.Set2(point.x, point.y);

        auto transformation = this->h1Space.GetElementTransformation(element.id);
        transformation->SetIntPoint(&integrationPoint);
        this->_density.GetGradient(*transformation, result);

        return {result[0], result[1]};
    }

    mfem::GridFunction convertH1toL2(const mfem::GridFunction &function) {
        mfem::BilinearForm A(&h1Space);
        A.AddDomainIntegrator(new mfem::MassIntegrator);
        A.Assemble();
        A.Finalize();

        mfem::MixedBilinearForm B(&l2Space, &h1Space);
        B.AddDomainIntegrator(new mfem::MassIntegrator);
        B.Assemble();
        B.Finalize();

        mfem::LinearForm b(&h1Space);
        B.Mult(function, b);

        mfem::DSmoother smoother(A.SpMat());

        mfem::GridFunction result(&h1Space);
        result = 0;
        mfem::PCG(A, smoother, b, result);
        return result;
    }
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
        if (!newIntersection)
            return newIntersection;
        newIntersection->orientation.direction = getDirection(intersection, laserRay.getCriticalDensity());
        return newIntersection;
    }

private:
    const MeshFunction &density;
    const GradientCalculator &gradientCalculator;

    Vector getDirection(const Intersection &intersection, Density criticalDensity) {
        const double n1 = std::sqrt(1 - density.getValue(*intersection.previousElement) / criticalDensity.asDouble);
        const double n2 = std::sqrt(1 - density.getValue(*intersection.nextElement) / criticalDensity.asDouble);

        const auto gradient = gradientCalculator.getGradient(intersection);
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
    DiscreteLine side{};
    side.segmentCount = 100;
    side.length = 1;
    auto mfemMesh = constructRectangleMesh(side, side);

    //auto mfemMesh = std::make_unique<mfem::Mesh>("test_mesh.vtk", 1, 0);


    mfem::L2_FECollection l2FiniteElementCollection(0, 2);
    mfem::H1_FECollection h1FiniteElementCollection(1, 2);
    mfem::FiniteElementSpace l2FiniteElementSpace(mfemMesh.get(), &l2FiniteElementCollection);
    mfem::FiniteElementSpace h1FiniteElementSpace(mfemMesh.get(), &h1FiniteElementCollection);

    mfem::GridFunction densityGridFunction(&l2FiniteElementSpace);
    mfem::FunctionCoefficient densityFunctionCoefficient(density);
    densityGridFunction.ProjectCoefficient(densityFunctionCoefficient);

    mfem::GridFunction absorbedEnergyGridFunction(&l2FiniteElementSpace);
    absorbedEnergyGridFunction = 0;
    //End of MFEM boilerplate ------------------------------------------------------------------------------------------

    Mesh mesh(mfemMesh.get());
    MeshFunction densityMeshFunction(densityGridFunction, l2FiniteElementSpace);
    MeshFunction absorbedEnergyMeshFunction(absorbedEnergyGridFunction, l2FiniteElementSpace);

    //ConstantGradientCalculator gradientCalculator(Vector(12.8e20, 0));
    H1GradientCalculator gradientCalculator(l2FiniteElementSpace, h1FiniteElementSpace);
    gradientCalculator.updateDensity(densityGridFunction);

    auto start = std::chrono::high_resolution_clock::now();

    Laser laser(
            Length{1315e-7},
            [](const Point) { return Vector(1, 0.2); },
            Gaussian(0.1),
            Point(-0.1, 0.2),
            Point(-0.1, 0.1)
    );

    laser.generateRays(100);
    laser.generateIntersections(
            mesh, SnellsLaw(densityMeshFunction, gradientCalculator),
            DontStop());
            //StopAtCritical(densityMeshFunction));

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
    //sol_sock << "solution\n" << *mfemMesh << absorbedEnergyGridFunction << std::flush;
    sol_sock << "solution\n" << *mfemMesh << densityGridFunction << std::flush;
}