
# OCUSF

OCUSF (Open channel unsteady flow) is a simulation of flow through an open channel (EG river) 

This program converted from BASIC to 'C' by K.Goode from original model in  "Water and Wasterwater Engineering Hydraulics"

## Prerequistes
1. g++
2. make
3. cppunit-devel (to build nad run tests)
4. gdb (to debug code , tests)

## Build
```console
make (or make debug)
```
Tested using g++ (GCC) 8.3.1 20190223 (Red Hat 8.3.1-2)
Open Visual Studio Code at root directory to run and debug code from ide

## Run from Commandline
Both input file and output files are optional  
```console
./ocusf.exe OCUSF.dat results.json
```
Results generated are the same as in book


## Build Tests
```console
testing> make (or make debug)
```
Open Visual Studio Code at testing directory to run and debug code from ide

## Run Unit Tests from Commandline

```console
testing>./ocusf_tests.exe
```

