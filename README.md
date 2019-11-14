# trace

> An executable to run raytracing simulations

## Dependencies
None for now.

git

## Compiling
The project uses cmake to compile. It will download and setup all its dependencies automatically.
As a user just do the following:
```
mdir build
cd build/
cmake ../
make
```
This will create an executable called "trace" inside the build directory.

## Usage

To run the executable do this from project root:
```
cd run/
../build/trace
```

This will create a file called "ray.json". To show the result call from run directory:
```
python3 ../scripts/plot_result.py mesh.json ray.json
```
Note that for this to work you will need python 3 installed.