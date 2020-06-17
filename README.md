# Puma: a microservice of extruded visualization (VTK-m) for fusion tokamak simulation (WDM)
## Why extruded mesh visualization?
XGC is a edge-fusion tokamak simulation code developed by the Department of Energy. It outputs what it calls an "rz" mesh, which be reconstructed into a 3D mesh such that:

```
auto pt = vtkm::Vec3f(r, phi, z);
```
where ```phi``` is 
```
auto phi = whichPlane * (vtkm::TwoPi() / NumberOfPlanes);
```

A tokamak is donut shaped, but in the previous code it's represented in cylindrical form rather than donut form. Further, there are a certain number of planes that represent the tokamak, and phi is dependent on which plane which "rz" point is currently being converted from the "rz" space to the XYZ-space. If we wanted it to be a donut shape, then the pt would be located in the XYZ-space as:

```
auto tp = vtkm::Vec3f(r * vtkm::Cos(phi), r * vtkm::Sin(phi), z);
```

Previously, converting the RZ-space representation to an XYZ-space representation meant building the complete triangular mesh. Now, VTK-m can construct the extruded mesh on-the-fly, without having to construct the full representation. This reduces bandwidth and memory costs.

## Installation
The easiest way to install Puma dependencies is to use [Spack](www.spack.io). We'll start there. First install spack.
```
git clone https://github.com/m-kim/wdmapp-spack
```
### Note: this is not the standard spack. How come?
Well, it's a fork of a fork. In other words, wdmapp has it's own spack repository, and Puma needed some customizations as well. What kind of customizations? Well, to keep VTK-m compiling smoothly, we opted to leave out the ```CELLSET_EXTRUDE``` from the ```CellSetList.``` in the ```Policy.``` For now, that means there needs to be a separate Puma repository at:
```
git clone https://gitlab.kitware.com/m-kim/vtk-m
``` 

on branch ```extrude-truncate.``` Don't worry though, ```Policy``` is getting a clean up after version 1.6, which should remove this obstacle. But, for now this is where things are. However, the custom spack fork (of a fork) will take care of that. 

### Spack is installed. Now what?
A few things need to be installed. First, activate spack:

```
. ~/wdm-spack/share/spack.setup.sh
```
Next, effis (and more importantly kittie) need to be installed. 

```
spack install effis@develop
```

There are a couple things that I should note. One, OpenMPI seems to be having problems with CMake right now. So, I usually like to install using ```^mpich``` but that's up to user discretion. Second, *running effis requires codar-cheetah.* It is not listed as running requirement, so will need to be installed by hand.

```
spack install codar-cheetah
```

There's a problem with that as well, because codar-cheetah needs python3. This can be problematic if you have code that needs python2 or running on rhea because it doesn't have python3. So, I (and [Eric](https://github.com/suchyta1/)) like to install local version of both pythons. I'll leave that as an exercise to the reader. 

Once those are installed, install the customized VTK-m:
```
spack install mvtkm
```

## Puma Installation
Now that all the pre-requisites are setup, let's install Puma:
```
git clone git@github.com:m-kim/Puma.git
```
and build it
```
spack load codar-cheetah
spack load effis
spack load mvtkm
cd Puma
mkdir build
cd build
cmake ../
make -j
```

## Running Puma
Congratulations, you've installed Puma. How do you run it? That's a great question, there are two example files in the ```effis-examples``` directory to help you out. 

### example.yaml
The first one is called, simply enough, ```example.yaml.``` I used this one to test at home using data generated on sdg-tm76. To use effis, first the ```example.yaml``` needs to be "composed." This sets up the directories and cheetah to run. 
```
effis-compose.py example.yaml
```

This will create directories in the *cwd* ```runs/run-1``` and pre-populate scripts for executing through cheetah. To run, the composed example needs to be "submitted."
```
effis-submit runs/run-1
```

Note, the directory that was created in the ```compose``` stage is now passed as *cli* to ```effis-submit.``` This will launch Puma in the background.

### rhea-example.yaml

Another effis composition file is ```rhea-example.yaml.``` This will launch Puma as a batch service through Cheetah (via slurm) on rhea, automatcially. The time out is set to ```3600```, and the charge code is setup to the SDG charge code. This looks for the specific files in a particular location on the fileserve.

#### Note about rhea
It's easier to install your own local pythons, both 2 and 3, and setup $SPACK_LOCATION/etc/packages.yaml to point to those directories. Also, CMake currently seems to have some issues with MPI (for example, [discourse discussion](https://discourse.cmake.org/t/cmake-cannot-find-mpi-in-standard-installation-path/850), [#18196](https://gitlab.kitware.com/cmake/cmake/-/issues/18196) ). My build workflow has been: install through spack ```effis@develop^mpich```, ```codar-cheetah``` and  ```mvtkm+mpi^mpich```. Then, I build ```puma```, which usually produces a MPI_C_FOUND, MPI_CXX_FOUND error). I run ```module load openmpi``` and ```module unload openmpi```, run cmake again and everything runs. ¯\\_(ツ)_/¯ This has been the most consistent approach to getting up and running on rhea.

### Puma specific yaml
Most of this should be explained in the effis documentation. Specifically for Puma, Puma needs to file names and path names. These should be passed in like this:


```
executable_path: /home/mark/Projects/Puma/build/wsl/Debug/pumacli/pumacli
    commandline_args: 
    - -meshname
    - xgc.mesh
    - -meshpathname
    - /mnt/c/Users/mark/Desktop/xgc-gene/
    - -diagname
    - xgc.oneddiag
    - -diagpathname
    - /mnt/c/Users/mark/Desktop/xgc-gene/
    - -filename
    - xgc.3d.00015
    - -filepathname
    - /mnt/c/Users/mark/Desktop/xgc-gene/
```
Note, the ```-meshname,``` ```-diagname``` and ```-filename``` variables do not end in ```.bp.``` 


## Output
Currently, the output is ```output-#.pnm``` i.e. pnm files. 

## Can I run Puma by itself?
No. Although Puma takes all it's files and paths in through the *cli*, it depends on kittie infrastructure to run. Kittie sets up a lot of bash variables in the background to ensure the workflow is composed correctly. It cannot be run by itself.

## No really, can I run Puma myself?
Yes, by going into the code and switching out the kittie for ADIOS calls.

# Conclusion
Puma is a visualization microservice for WDM using effis/kittie. By integrating kittie directly into Puma, it should allow for more flexible workflow management.
