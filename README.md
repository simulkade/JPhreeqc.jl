# Julia wrapper for PhreeqcRM package
This package wraps all the C functions of the PhreeqcRM package, except the MPI
functions (I did not know how to do it).

## Pease cite as
[![DOI](https://zenodo.org/badge/58193697.svg)](https://zenodo.org/badge/latestdoi/58193697)

## Installation
Open Julia and type
```
]add https://github.com/simulkade/JPhreeqc.jl
build JPhreeqc
```
The `build` command will download the latest binaries on Windows and Linux. Note that on windows, a compiled dll is provided in the `lib` folder. In addition, you need to have `Visual Studio 2019` with C++ tools installed for all the library dependencies. You can download [the community edition](https://www.visualstudio.com/en-us/products/visual-studio-community-vs.aspx) for free.  
I do not have a Mac so you are on your own there. Please feel free to create a pull request if you intend to make it work on Mac.

## Plan
  - ~~add a test functions~~
  - ~~connect and use it with JFVM.jl package~~
  - ~~Make the syntax more convenient, particularly regarding the usage of Int and Int32~~
  - Wrap IPhreeqc functions as well
  - ~~Write some convenience functions ~~
    * a `phreeqcrm` type?
  - add more tests for cases without transport
  - Convert some of the phreeqc examples to JPhreeqc

## About
This package is written by Ali A. Eftekhari. I use my personal time for its development, but I use it in my work at DHRTC.  
[PhreeqcRM](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and other Phreeqc programs are developed and distributed by U.S. Geological Survey (USGS).
