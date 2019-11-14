# trace

> An executable to run raytracing simulations

## Dependencies
To be able to compile and run this project you will need **git** (even <u>during compilation</u>
as it might be used to fetch some dependencies). Beside that you must have **cmake** and
a **c++ compiler** of your choice that supports **C++14**. The project is mainly developed
and tested using **g++**.

To ensure that dependencies are met you can do something along these lines depending on your distribution
(example is for ubuntu like distributions):
```
sudo apt-get update && sudo apt-get install git build-essential cmake
```
Then you will need the actual source which can be obtained simply by doing:
```
git clone https://github.com/SachCZ/trace
```


## Compiling
The project uses cmake to compile. It will download and setup all its dependencies automatically.
First go to the project root. If you followed the commands from previous paragraph just do this:
```
cd trace/
```
Then to actually compile the project execute these commands:
```
mkdir build
cd build/
cmake ../
make
```
This will create an executable called "trace" inside the build directory.

## Usage

To run the executable do the following from the project root:
```
cd run/
../build/trace
```

Files called "ray.json" and "mesh.json" will be generated. Up to this point the whole process
of compiling and using the code was tested on Ubuntu Server 18.04 LTS instance.

To show the result call this python script from the run directory:
```
python3 ../scripts/plot_result.py mesh.json ray.json
```
Note that for this to work you will need python 3 with matplotlib and numpy installed.

Now to actually do something useful please feel free to fiddle around with *main.cpp* in the
root directory as it is the entry point to the application. It contains 
a simple demonstration of the [raytracer library](https://github.com/SachCZ/raytracer)
to get you started.

