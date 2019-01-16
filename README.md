MadingleyCPP-OpenMP
Parallel (OpenMP) C++ version of the Madingley model, includes: 
- Minor bug fixes
- Functions for cohort consumption outputs
- Performance gain (30%), by reducing the number of prey considered by predators  

The computation time in the parallel version of the model scales with the number of threads as followed:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{time(s)}{timestep}=1000*exp(-0.032*N_{threads})"/>

For a global simulation at 1 degree and a maximum of 1000 cohorts in a single grid cell.

# Madingley C++ Guide (Linux/Ubuntu) 
A brief tutorial for getting the parallel C++ version of Madingley up and running on a linux machine

### 1. Prerequisites
- compiler with OpenMP support
- netCDF 4.6.0+
- netCDF-cxx4 4.2.1+

On a Ubuntu 16.04 machine the netCDF libraries can be installed using:
```bash
apt-get install libnetcdf-dev libnetcdf-cxx-legacy-dev libnetcdf-c++4-dev
```
The following commands can be used to check if both libraries are installed correctly:
```bash
nc-config --all
ncxx4-config --all
```
### 2. Clone the MadingleyCPP-opemMP repository
Currently the repository is set to private:
```bash
git clone https://github.com/SHoeks/MadingleyCPP-opemMP 
```
The **src** directory contains the source code, where the various functions required to run Madingley are distributed across 5 subdirectories: 
- Input, functions for reading model inputs
- Data, function for correctly setting spatial and temporal data
- Model, functions mostly relevant to the model ecology 
- Output, functions related to writing simulation outputs 
- Tools, functions for calculating variable values and sampling values from various distributions

Additionally, the **src** directory contains the folder **Model_setup**, which contains the setup .csv files (see **4. Running MadingleyCPP-opemMP**).  

### 3. Compile MadingleyCPP-ll
In addition to the source code, the src directory contains the **Makefile**, which expects all libraries to be installed in their default locations. If for some reason the libraries are located in another location, line 8 and 11 of the **Makefile** must include the full path of the corrosponding library. The 2 config commands shown under **Prerequisites** can be used to find the correct path of both netCDF libraries. If specified correctly the **Makefile** can be executed from the **src** directory: 
```bash
cd /MadingleyCPP-opemMP/src/
make
```

### 4. Running MadingleyCPP-opemMP
During compilation of the source code two new folders will be created: 
- /MadingleyCPP-opemMP/dist/
- /MadingleyCPP-opemMP/build/

Besides the Madingley executable the **dist** folder will contain the subdirectories:
- /MadingleyCPP-opemMP/dist/input
- /MadingleyCPP-opemMP/dist/output 

The **/MadingleyCPP-opemMP/dist/input** folder contains the following simulation setup files:
- SimulationControlParameters.csv, used to set simulation properties
- CohortFunctionalGroupDefinitions.csv, used to set initialization properties cohorts
- StockFunctionalGroupDefinitions.csv, used to set initialization properties stocks
- EnvironmentalDataLayers.csv, used to define the environmental data layers 
- MassBinDefinitions.csv, definition of the body mass bins
- OutputControlParameters.csv, used to select simulation outputs

In the SimulationControlParameters.csv the path the environmetal data (netCDF files) must be set (line 2). The environmental data can be downloaded from: https://link.to.netcdf/files


The **/MadingleyCPP-opemMP/dist/output** will be used to write simulation results to.

The **dist** folder contains the executable **madingley**, which can be run using:
```bash
cd /MadingleyCPP-opemMP/dist/
./madingley
```

# Madingley C++ Docker image (MacOS)
The simplify to process of getting the various dependencies installed correctly the following Docker image could be deployed: https://github.com/SHoeks/MadingleyCPP_Dockerfile

See https://docs.docker.com/docker-for-mac/ for a brief tutorial on installing Docker on a Mac.

After installing the Docker desktop application the Docker image can be downloaded using:
```bash
mkdir DockerMadingley
cd DockerMadingley
git clone https://SHoeks@github.com/SHoeks/MadingleyCPP_Dockerfile
```
The next step would be to build the image
```bash
docker build MadingleyCPP_Dockerfile
```
note the id after a successful build, for example: 
```bash
Successfully built e26ad9346bcd
```
The following command can be used to check if the image is now available:
```bash
docker images
```
The setupDocker bash script can now be used to setup the image and copy all required data. 
```bash
cd MadingleyCPP_Dockerfile
sh setupDocker.sh
```
The bash script will ask for the location of the environmental data, which can be downloaded from: https://link.to.netcdf/files. An example path the to folder would be:
```bash
/Users/user/Projects/MadingleyData-master/
```
In addition it will ask the name of the docker image, for example (from the example above):
```bash
e26ad9346bcd
```
Next, the setupDocker bash script will start the Docker image and mount the data directory. The terminal will now show an Ubuntu 16.04 shell. 

Check:
```bash
ls /mnt/MadingleyData-master/
```
to see if the data folder is mounted correctly.
```bash
nc-config --all
ncxx4-config --all
```
to see if the Netcdf libraries are installed correctly.
```bash
lscpu
```
to see if the hardware is recognized correctly.

The getMadingley bash script (which can be found from within the Docker image) can be used to get and compile the latest version from github. 
```bash
sh /home/getMadingley.sh
```
The compiled model can be found in the home folder.

# MadingleyCPP-opemMP Guide (Windows 10)
For Windows *'bash for Windows'* is available directly from Microsoft: https://docs.microsoft.com/en-us/windows/wsl/install-win10

The same steps as described for Linux/Ubuntu or MacOS can be applied from the Ubuntu bash shell within Windows.
