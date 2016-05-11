# Julia wrapper for PhreeqcRM package
This package wraps all the C functions of the PhreeqcRM package, except the MPI
functions (I did not know how to do it).  

## Installation
Open Julia and type `Pkg.clone(https://github.com/simulkade/JPhreeqc.jl.git)`

## Plan
  - ~~add a test functions~~
  - ~~connect and use it with JFVM.jl package~~
  - Make the syntax more convenient, particularly regarding the usage of Int and Int32
  - Wrap IPhreeqc functions as well
  - Write some convenience functions, and perhaps a `phreeqcrm` type
  - add more tests for cases without transport

## About
This package is written by Ali Akbar Eftekhari. I use my personal time for its development, but I use it in my work at DHRTC.  
[PhreeqcRM](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and other Phreeqc programs are developed and distributed by U.S. Geological Survey (USGS).
