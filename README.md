# gsIncompressibleFlow

An optional sumbodule of G+Smo providing tools for solving incompressible flow problems.

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsIncompressibleFlow"```|
|--:|---|
|License|![GitHub License](https://img.shields.io/github/license/gismo/gismo?color=008A00)|
|OS support|Linux, Windows, macOS|
|Developers/maintainers| [![Static Badge](https://img.shields.io/badge/@hhornik-008A00)](https://github.com/hhornik) [![Static Badge](https://img.shields.io/badge/@turnerov-008A00)](https://github.com/turnerov)

<video width="500" controls>
  <source src="https://raw.githubusercontent.com/gismo/gsIncompressibleFlow/media/video/profile_vel_unst.mp4" type="video/mp4">
</video>

## Installation

#### 1. Download G+smo

Similarly to other submodules, it is necessary to download G+Smo before downloading gsIncompressibleFlow:
```
git clone https://github.com/gismo/gismo.git
```

#### 2. Download gsIncompressibleFlow

To get the gsIncompressibleFlow submodule, configure G+Smo with the option `-DGISMO_OPTIONAL="<other submodules>;gsIncompressibleFlow"` (in addition to any other CMake options):
```
cd gismo
mkdir build
cd build
cmake .. -DGISMO_OPTIONAL="<other submodules>;gsIncompressibleFlow"
```
This will trigger a download of gsIncompressibleFlow from GitHub. 

Since optional submodules of G+Smo are independent .git repositories, you need to access them via their own repositories located in `/path/to/gismo/optional/` to get different versions/branches or keep your copy of the submodule up-to-date. For example, to obtain the latest version of gsIncompressibleFlow, do
```
cd /path/to/gismo/optional/gsIncompressibleFlow
git pull
```

#### 3. Compile

Once gsIncompressibleFlow is downloaded, it is compiled together with G+Smo:
```
cd /path/to/gismo/build
make
```
When the compilation is complete, you can find the compiled library in `/path/to/gismo/build/lib` and examples in `/path/to/gismo/build/bin`.


## Usage

The usage and capabilities of the submodule are demonstrated by examples in `/path/to/gismo/optional/gsIncompressibleFlow/examples`. To display the usage information, do (for example):
```
cd /path/to/gismo/build/bin
./gsINSSolversExample -h
```
