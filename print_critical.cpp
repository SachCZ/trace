#include <iostream>
#include <raytracer.h>

int main(int argc, char* argv[]){
    auto wavelength = raytracer::Length{std::stod(argv[1])};
    std::cout << raytracer::calcCritDens(wavelength).asDouble;
}

