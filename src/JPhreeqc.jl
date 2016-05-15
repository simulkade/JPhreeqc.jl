module JPhreeqc

if OS_NAME==:Windows
  const Lib_PhreeqcRM_path=joinpath(Pkg.dir(), "JPhreeqc", "deps", "PhreeqcRMd.dll")
else
  const Lib_PhreeqcRM_path=joinpath("/usr/local/lib", "libphreeqcrm.so")
end
#const phreeqc_path1=Libdl.find_library(["libphreeqcrm"], ["/usr/local/lib"])

# Note: int *a -> type: Ptr{Cint} value: ones(Cint, 5)
# Note: double *c -> type: Ptr{Cdouble} value ones(Cdouble, 5)
# Note: char* (string) -> Cstring or Ptr{UInt8}
"""
Creates a reaction module. If the code is compiled with the preprocessor directive
USE_OPENMP, the reaction module is multithreaded. If the code is compiled
with the preprocessor directive USE_MPI, the reaction module will use MPI
and multiple processes. If neither preprocessor directive is used,
the reaction module will be serial (unparallelized).

Parameters
nxyz	The number of grid cells in the user's model.
n_threads	(or comm, MPI) When using OPENMP, the argument (nthreads) is
the number of worker threads to be used. If nthreads <= 0, the number of
threads is set equal to the number of processors of the computer. When using MPI,
the argument (comm) is the MPI communicator to use within the reaction module.
Return values
Id	of the PhreeqcRM instance, negative is failure (See RM_DecodeError).
See also
RM_Destroy
"""
function RM_Create(nxyz::Int, n_threads::Int)
  phreeqc_id= ccall((:RM_Create, Lib_PhreeqcRM_path), Cint, (Cint, Cint), nxyz, n_threads)
  return convert(Int, phreeqc_id)
end

"""
Close the output and log files.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_OpenFiles, RM_SetFilePrefix
C Example:
 status = RM_CloseFiles(id);
MPI:
Called only by root.
"""
function RM_CloseFiles(id::Int)
  IRM_RESULT=ccall((:RM_CloseFiles, Lib_PhreeqcRM_path), Cint, (Cint,), id)
end

"""
N sets of component concentrations are converted to SOLUTIONs numbered 1-n in the Utility IPhreeqc. The solutions can be reacted and manipulated with the methods of IPhreeqc. If solution concentration units (RM_SetUnitsSolution) are per liter, one liter of solution is created in the Utility instance; if solution concentration units are mass fraction, one kilogram of solution is created in the Utility instance. The motivation for this method is the mixing of solutions in wells, where it may be necessary to calculate solution properties (pH for example) or react the mixture to form scale minerals. The code fragments below make a mixture of concentrations and then calculate the pH of the mixture.

Parameters
id	The instance id returned from RM_Create.
c	Array of concentrations to be made SOLUTIONs in Utility IPhreeqc. Array storage is equivalent to Fortran (n,ncomps).
n	The number of sets of concentrations.
tc	Array of temperatures to apply to the SOLUTIONs, in degree C. Array of size n.
p_atm	Array of pressures to apply to the SOLUTIONs, in atm. Array of size n.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 c_well = (double*) malloc((size_t) ((size_t) (1 * ncomps * sizeof(double))));
for (i = 0; i < ncomps; i++)
{
  c_well[i] = 0.5 * c[0 + nxyz*i] + 0.5 * c[9 + nxyz*i];
}
tc = (double*) malloc((size_t) (1 * sizeof(double)));
p_atm = (double*) malloc((size_t) (1 * sizeof(double)));
tc[0] = 15.0;
p_atm[0] = 3.0;
iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm);
strcpy(str, "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1");
status = RunString(iphreeqc_id, str);
status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5);
status = GetSelectedOutputValue2(iphreeqc_id, 1, 0, &vtype, &pH, svalue, 100);
"""
function RM_Concentrations2Utility(id::Int, c::Array{Float64,1}, n::Int,
  tc::Array{Float64,1}, p_atm::Array{Float64,1})
  IRM_RESULT=ccall((:RM_Concentrations2Utility, Lib_PhreeqcRM_path),
  Cint, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), id, c, n, tc, p_atm)
end

"""
Provides a mapping from grid cells in the user's model to reaction cells in PhreeqcRM. The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells for which chemistry must be run. The array of size nxyz (the number of grid cells, RM_GetGridCellCount) must contain the set of all integers 0 <= i < count_chemistry, where count_chemistry is a number less than or equal to nxyz. Inactive cells are assigned a negative integer. The mapping may be many-to-one to account for symmetry. Default is a one-to-one mapping–all user grid cells are reaction cells (equivalent to grid2chem values of 0,1,2,3,...,nxyz-1).

Parameters
id	The instance id returned from RM_Create.
grid2chem	An array of integers: Nonnegative is a reaction cell number (0 based), negative is an inactive cell. Array of size nxyz (number of grid cells).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 // For demonstation, two equivalent rows by symmetry
grid2chem = (int*) malloc((size_t) (nxyz * sizeof(int)));
for (i = 0; i < nxyz/2; i++)
{
  grid2chem[i] = i;
  grid2chem[i+nxyz/2] = i;
}
status = RM_CreateMapping(id, grid2chem);
"""
function RM_CreateMapping(id::Int, grid2chem::Array{Int,1})
  grid2chem=convert(Array{Int32}, grid2chem)
  IRM_RESULT=ccall((:RM_CreateMapping, Lib_PhreeqcRM_path),
  Cint, (Cint, Ptr{Cint}), id, grid2chem)
end

"""
Abort the program. Result will be interpreted as an IRM_RESULT value and decoded;
err_str will be printed; and the reaction module will be destroyed. If using MPI, an MPI_Abort message will be sent before the reaction module is destroyed. If the id is an invalid instance, RM_Abort will return a value of IRM_BADINSTANCE, otherwise the program will exit with a return code of 4.

Parameters
id	The instance id returned from RM_Create.
result	Integer treated as an IRM_RESULT return code.
err_str	String to be printed as an error message.
Return values
IRM_RESULT	Program will exit before returning unless id is an invalid reaction module id.
See also
RM_Destroy, RM_ErrorMessage.
C Example:
 iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm);
strcpy(str, "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1");
status = RunString(iphreeqc_id, str);
if (status != 0) status = RM_Abort(id, status, "IPhreeqc RunString failed");
"""
function RM_Abort(id::Int, result::Int, err_str::AbstractString)
  IRM_RESULT=ccall((:RM_Abort, Lib_PhreeqcRM_path), Cint,
  (Cint, Cint, Cstring), id, result, err_str)
end

"""
Destroys a reaction module.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_Create
C Example:
 status = RM_Destroy(id);
"""
function RM_Destroy(id::Int)
  IRM_RESULT=ccall((:RM_Destroy, Lib_PhreeqcRM_path), Cint,
  (Cint,), id)
end

"""
If e is negative, this method prints an error message corresponding to IRM_RESULT e. If e is non-negative, no action is taken.

Parameters
id	The instance id returned from RM_Create.
e	An IRM_RESULT value returned by one of the reaction-module methods.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
IRM_RESULT definition:
 typedef enum {
  IRM_OK            =  0,  //Success
  IRM_OUTOFMEMORY   = -1,  //Failure, Out of memory
  IRM_BADVARTYPE    = -2,  //Failure, Invalid VAR type
  IRM_INVALIDARG    = -3,  //Failure, Invalid argument
  IRM_INVALIDROW    = -4,  //Failure, Invalid row
  IRM_INVALIDCOL    = -5,  //Failure, Invalid column
  IRM_BADINSTANCE   = -6,  //Failure, Invalid rm instance id
  IRM_FAIL          = -7,  //Failure, Unspecified
} IRM_RESULT;
C Example:
 status = RM_CreateMapping(id, grid2chem);
if (status < 0) status = RM_DecodeError(id, status);
"""
function RM_DecodeError(id::Int, e_result::Int)
  IRM_RESULT=ccall((:RM_DecodeError,Lib_PhreeqcRM_path),Cint,
  (Cint,Cint),id,e_result)
end

"""
Send an error message to the screen, the output file, and the log file.

Parameters
id	The instance id returned from RM_Create.
errstr	String to be printed.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_OpenFiles, RM_LogMessage, RM_OutputMessage, RM_ScreenMessage, RM_WarningMessage.
C Example:
 status = RM_ErrorMessage(id, "Goodbye world");
"""
function RM_ErrorMessage(id::Int, errstr::AbstractString)
  IRM_RESULT=ccall((:RM_ErrorMessage,Lib_PhreeqcRM_path),Cint,(Cint,Cstring),id,errstr)
end

"""
Returns the number of items in the list of all elements in the InitialPhreeqc instance. Elements are those that have been defined in a solution or any other reactant (EQUILIBRIUM_PHASE, KINETICS, and others). The method can be called multiple times and the list that is created is cummulative. The list is the set of components that needs to be transported. By default the list includes water, excess H and excess O (the H and O not contained in water); alternatively, the list may be set to contain total H and total O (RM_SetComponentH2O), which requires transport results to be accurate to eight or nine significant digits. If multicomponent diffusion (MCD) is to be modeled, there is a capability to retrieve aqueous species concentrations (RM_GetSpeciesConcentrations) and to set new solution concentrations after MCD by using individual species concentrations (RM_SpeciesConcentrations2Module). To use these methods the save-species property needs to be turned on (RM_SetSpeciesSaveOn). If the save-species property is on, RM_FindComponents will generate a list of aqueous species (RM_GetSpeciesCount, RM_GetSpeciesName), their diffusion coefficients at 25 C (RM_GetSpeciesD25), their charge (RM_GetSpeciesZ).

Parameters
id	The instance id returned from RM_Create.
Return values
Number	of components currently in the list, or IRM_RESULT error code (see RM_DecodeError).
See also
RM_GetComponent, RM_SetSpeciesSaveOn, RM_GetSpeciesConcentrations, RM_SpeciesConcentrations2Module, RM_GetSpeciesCount, RM_GetSpeciesName, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_SetComponentH2O.
C Example:
 // Get list of components
 ncomps = RM_FindComponents(id);
components = (char**) malloc((size_t) (ncomps* sizeof(char*)));
for (i = 0; i < ncomps; i++)
{
  components[i] = (char*) malloc((size_t) (100 * sizeof(char*)));
  status = RM_GetComponent(id, i, components[i], 100);
}
"""
function RM_FindComponents(id::Int)
  n_components=ccall((:RM_FindComponents, Lib_PhreeqcRM_path), Cint,
  (Cint,), id)
  return convert(Int, n_components)
end

"""
Fills an array with the cell numbers in the user's numbering sytstem that map to a cell in the PhreeqcRM numbering system. The mapping is defined by RM_CreateMapping.

Parameters
id	The instance id returned from RM_Create.
n	A cell number in the PhreeqcRM numbering system (0 <= n < RM_GetChemistryCellCount).
list	Array to store the user cell numbers mapped to PhreeqcRM cell n.
size	Input, the allocated size of list; it is an error if the array is too small. Output, the number of cells mapped to cell n.
Return values
IRM_RESULT	error code (see RM_DecodeError).
See also
RM_CreateMapping, RM_GetChemistryCellCount, RM_GetGridCellCount.
C Example:
 if (RM_GetBackwardMapping(rm_id, rm_cell_number, list, &size) == 0)
{
  if (strcmp(str, "HYDRAULIC_K") == 0)
  {
    return data->K_ptr[list[0]];
  }
}
"""
function RM_GetBackwardMapping(id::Int, n::Int, list::Array{Int32,1}, size1::Array{Int32,1})
  IRM_RESULT=ccall((:RM_GetBackwardMapping, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Ptr{Cint},Ptr{Cint}), id,n,list,size1)
end

"""
Returns the number of chemistry cells in the reaction module. The number of chemistry cells is defined by the set of non-negative integers in the mapping from user grid cells (RM_CreateMapping). The number of chemistry cells is less than or equal to the number of cells in the user's model.

Parameters
id	The instance id returned from RM_Create.
Return values
Number	of chemistry cells, or IRM_RESULT error code (see RM_DecodeError).
See also
RM_CreateMapping, RM_GetGridCellCount.
C Example:
 status = RM_CreateMapping(id, grid2chem);
nchem = RM_GetChemistryCellCount(id);
"""
function RM_GetChemistryCellCount(id::Int)
  n_chem_cells=ccall((:RM_GetChemistryCellCount, Lib_PhreeqcRM_path), Cint,
  (Cint,), id)
  return convert(Int, n_chem_cells)
end

"""
Retrieves an item from the reaction-module component list that was generated by calls to RM_FindComponents.

Parameters
id	The instance id returned from RM_Create.
num	The number of the component to be retrieved. C, 0 based.
chem_name	The string value associated with component num.
l	The length of the maximum number of characters for chem_name.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetComponentCount
C Example:
 // Get list of components
ncomps = RM_FindComponents(id);
components = (char**) malloc((size_t) (ncomps * sizeof(char*)));
for (i = 0; i < ncomps; i++)
{
  components[i] = (char*) malloc((size_t) (100* sizeof(char*)));
  status = RM_GetComponent(id, i, components[i], 100);
}
"""
function RM_GetComponent(id::Int, num::Int, chem_name::AbstractString, l::Int)
  IRM_RESULT=ccall((:RM_GetComponent, Lib_PhreeqcRM_path), Cint,
  (Cint, Cint, Ptr{UInt8}, Cint), id, num, chem_name, l)
end

"""
Returns the number of components in the reaction-module component list. The component list is generated by calls to RM_FindComponents. The return value from the last call to RM_FindComponents is equal to the return value from RM_GetComponentCount.

Parameters
id	The instance id returned from RM_Create.
Return values
The	number of components in the reaction-module component list, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetComponent.
C Example:
 ncomps1 = RM_GetComponentCount(id);
"""
function RM_GetComponentCount(id::Int)
  n_components=ccall((:RM_GetComponentCount, Lib_PhreeqcRM_path), Cint, (Cint,), id)
  return convert(Int, n_components)
end

"""
Transfer solution concentrations from each reaction cell to the concentration array given in the argument list (c). Units of concentration for c are defined by RM_SetUnitsSolution. For concentration units of per liter, the solution volume is used to calculate the concentrations for c. For mass fraction concentration units, the solution mass is used to calculate concentrations for c. Two options are available for the volume and mass of solution that are used in converting to transport concentrations: (1) the volume and mass of solution are calculated by PHREEQC, or (2) the volume of solution is the product of saturation (RM_SetSaturation), porosity (RM_SetPorosity), and representative volume (RM_SetRepresentativeVolume), and the mass of solution is volume times density as defined by RM_SetDensity. RM_UseSolutionDensityVolume determines which option is used. For option 1, the databases that have partial molar volume definitions needed to accurately calculate solution volume are phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
c	Array to receive the concentrations. Dimension of the array is equivalent to Fortran (nxyz, ncomps), where nxyz is the number of user grid cells and ncomps is the result of RM_FindComponents or RM_GetComponentCount. Values for inactive cells are set to 1e30.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetComponentCount, RM_GetSaturation, RM_SetConcentrations, RM_SetDensity, RM_SetRepresentativeVolume, RM_SetSaturation, RM_SetUnitsSolution, RM_UseSolutionDensityVolume.
C Example:
 c = (double*) malloc((size_t) (ncomps * nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);
"""
function RM_GetConcentrations(id::Int, c::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetConcentrations, Lib_PhreeqcRM_path), Cint,
  (Cint, Ptr{Cdouble}), id, c)
end

"""
Transfer solution densities from the reaction cells to the array given in the argument list (density). Densities are those calculated by the reaction module. Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate density: phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
density	Array to receive the densities. Dimension of the array is nxyz, where nxyz is the number of user grid cells (RM_GetGridCellCount). Values for inactive cells are set to 1e30.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 density = (double*) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetDensity(id, density);
"""
function RM_GetDensity(id::Int, density::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetDensity, Lib_PhreeqcRM_path), Cint,
  (Cint, Ptr{Cdouble}), id, density)
end

"""
Returns an array with the ending cell numbers from the range of cell numbers assigned to each worker.

Parameters
id	The instance id returned from RM_Create.
ec	Array to receive the ending cell numbers. Dimension of the array is the number of threads (OpenMP) or the number of processes (MPI).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_Create, RM_GetMpiTasks, RM_GetStartCell, RM_GetThreadCount.
C Example:
 n = RM_GetThreadCount(id) * RM_GetMpiTasks(id);
ec = (int*) malloc((size_t) (n * sizeof(int)));
status = RM_GetEndCell(id, ec);
"""
function RM_GetEndCell(id::Int, ec::Array{Int32,1})
  IRM_RESULT=ccall((:RM_GetEndCell, Lib_PhreeqcRM_path), Cint, (Cint, Ptr{Cint}), id, ec)
end

"""
Returns a string containing error messages related to the last call to a PhreeqcRM method to the character argument (errstr).

Parameters
id	The instance id returned from RM_Create.
errstr	The error string related to the last call to a PhreeqcRM method.
l	Maximum number of characters that can be written to the argument (errstr).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 if (status != IRM_OK)
{
  l = RM_GetErrorStringLength(id);
  errstr = (char*) malloc((size_t) (l * sizeof(char) + 1));
  RM_GetErrorString(id, errstr, l+1);
  fprintf(stderr,"%s", errstr);
  free(errstr);
  RM_Destroy(id);
  exit(1);
}
"""
function RM_GetErrorString(id::Int, errstr::AbstractString, l::Int)
  IRM_RESULT=ccall((:RM_GetErrorString, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring,Cint), id, errstr, l)
end

"""
Returns the length of the string that contains error messages related to the last call to a PhreeqcRM method.

Parameters
id	The instance id returned from RM_Create.
Return values
int	Length of the error message string (for C, equivalent to strlen, does not include terminating \0).
C Example:
 if (status != IRM_OK)
{
  l = RM_GetErrorStringLength(id);
  errstr = (char*) malloc((size_t) (l * sizeof(char) + 1));
  RM_GetErrorString(id, errstr, l+1);
  fprintf(stderr,"%s", errstr);
  free(errstr);
  RM_Destroy(id);
  exit(1);
}
"""
function RM_GetErrorStringLength(id::Int)
  er_length=ccall((:RM_GetErrorStringLength, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the reaction-module file prefix to the character argument (prefix).

Parameters
id	The instance id returned from RM_Create.
prefix	Character string where the prefix is written.
l	Maximum number of characters that can be written to the argument (prefix).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetFilePrefix.
C Example:
 char str[100], str1[200];
status = RM_GetFilePrefix(id, str, 100);
strcpy(str1, "File prefix: ");
strcat(str1, str);
strcat(str1, "\n");
status = RM_OutputMessage(id, str1);
"""
function RM_GetFilePrefix(id::Int, prefix::AbstractString, l::Int)
  IRM_RESULT=ccall((:RM_GetFilePrefix, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring,Cint), id, prefix, l)
end

"""
Returns the gram formula weights (g/mol) for the components in the reaction-module component list.

Parameters
id	The instance id returned from RM_Create.
gfw	Array to receive the gram formula weights. Dimension of the array is ncomps, where ncomps is the number of components in the component list.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetComponentCount, RM_GetComponent.
C Example:
 ncomps = RM_FindComponents(id);
components = (char**) malloc((size_t) (ncomps * sizeof(char*)));
gfw = (double*) malloc((size_t) (ncomps * sizeof(double)));
status = RM_GetGfw(id, gfw);
for (i = 0; i < ncomps; i++)
{
  components[i] = (char*) malloc((size_t) (100 * sizeof(char*)));
  status = RM_GetComponent(id, i, components[i], 100);
  sprintf(str,"%10s    %10.3f\n", components[i], gfw[i]);
  status = RM_OutputMessage(id, str);
}
"""
function RM_GetGfw(id::Int, gfw::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetGfw, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, gfw)
end

"""
Returns the number of grid cells in the user's model, which is defined in the call to RM_Create. The mapping from grid cells to reaction cells is defined by RM_CreateMapping. The number of reaction cells may be less than the number of grid cells if there are inactive regions or symmetry in the model definition.

Parameters
id	The instance id returned from RM_Create.
Return values
Number	of grid cells in the user's model, negative is failure (See RM_DecodeError).
See also
RM_Create, RM_CreateMapping.
C Example:
 nxyz = RM_GetGridCellCount(id);
sprintf(str1, "Number of grid cells in the user's model: %d\n", nxyz);
status = RM_OutputMessage(id, str1);
"""
function RM_GetGridCellCount(id::Int)
  n_cells=ccall((:RM_GetGridCellCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns an IPhreeqc id for the ith IPhreeqc instance in the reaction module.

For the threaded version, there are nthreads + 2 IPhreeqc instances, where nthreads is defined in the constructor (RM_Create). The number of threads can be determined by RM_GetThreadCount. The first nthreads (0 based) instances will be the workers, the next (nthreads) is the InitialPhreeqc instance, and the next (nthreads + 1) is the Utility instance. Getting the IPhreeqc pointer for one of these instances allows the user to use any of the IPhreeqc methods on that instance. For MPI, each process has exactly three IPhreeqc instances, one worker (number 0), one InitialPhreeqc instance (number 1), and one Utility instance (number 2).

Parameters
id	The instance id returned from RM_Create.
i	The number of the IPhreeqc instance to be retrieved (0 based).
Return values
IPhreeqc	id for the ith IPhreeqc instance, negative is failure (See RM_DecodeError).
See also
RM_Create, RM_GetThreadCount. See IPhreeqc documentation for descriptions of IPhreeqc methods.
C Example:
 // Utility pointer is worker number nthreads + 1
iphreeqc_id1 = RM_GetIPhreeqcId(id, RM_GetThreadCount(id) + 1);
"""
function RM_GetIPhreeqcId(id::Int, i::Int)
  iphreeqc_id=ccall((:RM_GetIPhreeqcId, Lib_PhreeqcRM_path), Cint,  (Cint,Cint), id, i)
end

"""
Returns the MPI task number. For the OPENMP version, the task number is always zero and the result of RM_GetMpiTasks is one. For the MPI version, the root task number is zero, and all workers have a task number greater than zero. The number of tasks can be obtained with RM_GetMpiTasks. The number of tasks and computer hosts are determined at run time by the mpiexec command, and the number of reaction-module processes is defined by the communicator used in constructing the reaction modules (RM_Create).

Parameters
id	The instance id returned from RM_Create.
Return values
The	MPI task number for a process, negative is failure (See RM_DecodeError).
See also
RM_GetMpiTasks.
C Example:
 sprintf(str1, "MPI task number:  %d\n", RM_GetMpiMyself(id));
status = RM_OutputMessage(id, str1);
"""
function RM_GetMpiMyself(id::Int)
  mpi_task_number=ccall((:RM_GetMpiMyself, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the number of MPI processes (tasks) assigned to the reaction module. For the OPENMP version, the number of tasks is always one (although there may be multiple threads, RM_GetThreadCount), and the task number returned by RM_GetMpiMyself is zero. For the MPI version, the number of tasks and computer hosts are determined at run time by the mpiexec command. An MPI communicator is used in constructing reaction modules for MPI. The communicator may define a subset of the total number of MPI processes. The root task number is zero, and all workers have a task number greater than zero.

Parameters
id	The instance id returned from RM_Create.
Return values
The	number of MPI processes assigned to the reaction module, negative is failure (See RM_DecodeError).
See also
RM_GetMpiMyself, RM_Create.
C Example:
 mpi_tasks = RM_GetMpiTasks(id);
sprintf(str1, "Number of MPI processes: %d\n", mpi_tasks);
status = RM_OutputMessage(id, str1);
"""
function RM_GetMpiTasks(id::Int)
  n_mpi_process=ccall((:RM_GetMpiTasks, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the user number for the nth selected-output definition. Definitions are sorted by user number. Phreeqc allows multiple selected-output definitions, each of which is assigned a nonnegative integer identifier by the user. The number of definitions can be obtained by RM_GetSelectedOutputCount. To cycle through all of the definitions, RM_GetNthSelectedOutputUserNumber can be used to identify the user number for each selected-output definition in sequence. RM_SetCurrentSelectedOutputUserNumber is then used to select that user number for selected-output processing.

Parameters
id	The instance id returned from RM_Create.
n	The sequence number of the selected-output definition for which the user number will be returned. C, 0 based.
Return values
The	user number of the nth selected-output definition, negative is failure (See RM_DecodeError).
See also
RM_GetSelectedOutput, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputCount, RM_GetSelectedOutputHeading, RM_GetSelectedOutputRowCount, RM_SetCurrentSelectedOutputUserNumber, RM_SetSelectedOutputOn.
C Example:
 for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  fprintf(stderr, "Selected output sequence number: %d\n", isel);
  fprintf(stderr, "Selected output user number:     %d\n", n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double*) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
"""
function RM_GetNthSelectedOutputUserNumber(id::Int, n::Int)
  nth_output=ccall((:RM_GetNthSelectedOutputUserNumber, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, n)
end

"""
Returns a vector of saturations (sat) as calculated by the reaction module. Reactions will change the volume of solution in a cell. The transport code must decide whether to ignore or account for this change in solution volume due to reactions. Following reactions, the cell saturation is calculated as solution volume (RM_GetSolutionVolume) divided by the product of representative volume (RM_SetRepresentativeVolume) and the porosity (RM_SetPorosity). The cell saturation returned by RM_GetSaturation may be less than or greater than the saturation set by the transport code (RM_SetSaturation), and may be greater than or less than 1.0, even in fully saturated simulations. Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate solution volume and saturation: phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
sat_calc	Vector to receive the saturations. Dimension of the array is set to nxyz, where nxyz is the number of user grid cells (RM_GetGridCellCount). Values for inactive cells are set to 1e30.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetSolutionVolume, RM_SetPorosity, RM_SetRepresentativeVolume, RM_SetSaturation.
C Example:
 sat_calc = (double*) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSaturation(id, sat_calc);
"""
function RM_GetSaturation(id::Int, sat_calc::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSaturation, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, sat_calc)
end

"""
Populates an array with values from the current selected-output definition. RM_SetCurrentSelectedOutputUserNumber determines which of the selected-output definitions is used to populate the array.

Parameters
id	The instance id returned from RM_Create.
so	An array to contain the selected-output values. Size of the array is equivalent to Fortran (nxyz, col), where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount), and col is the number of columns in the selected-output definition (RM_GetSelectedOutputColumnCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetNthSelectedOutputUserNumber, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputCount, RM_GetSelectedOutputHeading, RM_GetSelectedOutputRowCount, RM_SetCurrentSelectedOutputUserNumber, RM_SetSelectedOutputOn.
C Example:
 for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double*) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
"""
function RM_GetSelectedOutput(id::Int, so::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSelectedOutput, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, so)
end

"""

"""
function RM_GetSelectedOutputColumnCount(id::Int)
  n_columns_output=ccall((:RM_GetSelectedOutputColumnCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the number of selected-output definitions. RM_SetCurrentSelectedOutputUserNumber determines which of the selected-output definitions is used.

Parameters
id	The instance id returned from RM_Create.
Return values
Number	of selected-output definitions, negative is failure (See RM_DecodeError).
See also
RM_GetNthSelectedOutputUserNumber, RM_GetSelectedOutput, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputHeading, RM_GetSelectedOutputRowCount, RM_SetCurrentSelectedOutputUserNumber, RM_SetSelectedOutputOn.
C Example:
 for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double*) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
"""
function RM_GetSelectedOutputCount(id::Int)
  n_so_count=ccall((:RM_GetSelectedOutputCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns a selected output heading. The number of headings is determined by RM_GetSelectedOutputColumnCount. RM_SetCurrentSelectedOutputUserNumber determines which of the selected-output definitions is used.

Parameters
id	The instance id returned from RM_Create.
icol	The sequence number of the heading to be retrieved. C, 0 based.
heading	A string buffer to receive the heading.
length	The maximum number of characters that can be written to the string buffer.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetNthSelectedOutputUserNumber, RM_GetSelectedOutput, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputCount, RM_GetSelectedOutputRowCount, RM_SetCurrentSelectedOutputUserNumber, RM_SetSelectedOutputOn.
C Example:
 char heading[100];
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  for (j = 0; j < col; j++)
  {
    status = RM_GetSelectedOutputHeading(id, j, heading, 100);
    fprintf(stderr, "          %2d %10s\n", j, heading);
  }
}
"""
function RM_GetSelectedOutputHeading(id::Int, icol::Int, heading::AbstractString, length::Int)
  IRM_RESULT=ccall((:RM_GetSelectedOutputHeading, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint, Cstring, Cint), id, icol, heading, length)
end

"""
Returns the number of rows in the current selected-output definition. However, the method is included only for convenience; the number of rows is always equal to the number of grid cells in the user's model, and is equal to RM_GetGridCellCount.

Parameters
id	The instance id returned from RM_Create.
Return values
Number	of rows in the current selected-output definition, negative is failure (See RM_DecodeError).
See also
RM_GetNthSelectedOutputUserNumber, RM_GetSelectedOutput, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputCount, RM_GetSelectedOutputHeading, RM_SetCurrentSelectedOutputUserNumber, RM_SetSelectedOutputOn.
C Example:
 for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double*) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Print results
  for (i = 0; i < RM_GetSelectedOutputRowCount(id)/2; i++)
  {
    fprintf(stderr, "Cell number %d\n", i);
    fprintf(stderr, "     Selected output: \n");
    for (j = 0; j < col; j++)
    {
      status = RM_GetSelectedOutputHeading(id, j, heading, 100);
      fprintf(stderr, "          %2d %10s: %10.4f\n", j, heading, selected_out[j*nxyz + i]);
    }
  }
  free(selected_out);
}
"""
function RM_GetSelectedOutputRowCount(id::Int)
  n_row_output=ccall((:RM_GetSelectedOutputRowCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Transfer solution volumes from the reaction cells to the array given in the argument list (vol). Solution volumes are those calculated by the reaction module. Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
vol	Array to receive the solution volumes. Dimension of the array is (nxyz), where nxyz is the number of user grid cells. Values for inactive cells are set to 1e30.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetSaturation.
C Example:
 volume = (double*) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSolutionVolume(id, volume);
"""
function RM_GetSolutionVolume(id::Int, vol::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSolutionVolume, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, vol)
end

"""
Transfer concentrations of aqueous species to the array argument (species_conc) This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. The list of aqueous species is determined by RM_FindComponents and includes all aqueous species that can be made from the set of components. Solution volumes used to calculate mol/L are calculated by the reaction module. Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
species_conc	Array to receive the aqueous species concentrations. Dimension of the array is (nxyz, nspecies), where nxyz is the number of user grid cells (RM_GetGridCellCount), and nspecies is the number of aqueous species (RM_GetSpeciesCount). Concentrations are moles per liter. Values for inactive cells are set to 1e30.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_GetSpeciesName, RM_SpeciesConcentrations2Module, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double*) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesConcentrations(id, species_c);
"""
function RM_GetSpeciesConcentrations(id::Int, species_conc::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSpeciesConcentrations, Lib_PhreeqcRM_path), Cint,  (
  Cint,Ptr{Cdouble}), id, species_conc)
end

"""
The number of aqueous species used in the reaction module. This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. The list of aqueous species is determined by RM_FindComponents and includes all aqueous species that can be made from the set of components.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	The number of aqueous species, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_GetSpeciesName, RM_SpeciesConcentrations2Module, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double*) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesConcentrations(id, species_c);
"""
function RM_GetSpeciesCount(id::Int)
  n_aqueous_species=ccall((:RM_GetSpeciesCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Transfers diffusion coefficients at 25C to the array argument (diffc). This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally in the database file. Databases distributed with the reaction module that have diffusion coefficients defined are phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
diffc	Array to receive the diffusion coefficients at 25 C, m^2/s. Dimension of the array is nspecies, where nspecies is is the number of aqueous species (RM_GetSpeciesCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesZ, RM_GetSpeciesName, RM_SpeciesConcentrations2Module, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
diffc = (double*) malloc((size_t) (nspecies * sizeof(double)));
status = RM_GetSpeciesD25(id, diffc);
"""
function RM_GetSpeciesD25(id::Int, diffc::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSpeciesD25, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, diffc)
end

"""
Transfers the name of the ith aqueous species to the character argument (name). This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. The list of aqueous species is determined by RM_FindComponents and includes all aqueous species that can be made from the set of components.

Parameters
id	The instance id returned from RM_Create.
i	Sequence number of the species in the species list. C, 0 based.
name	Character array to receive the species name.
length	Maximum length of string that can be stored in the character array.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_SpeciesConcentrations2Module, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 char name[100];
...
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
for (i = 0; i < nspecies; i++)
{
  status = RM_GetSpeciesName(id, i, name, 100);
  fprintf(stderr, "%s\n", name);
}
"""
function RM_GetSpeciesName(id::Int, i::Int, name::AbstractString, length::Int)
  IRM_RESULT=ccall((:RM_GetSpeciesName, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Cstring,Cint), id, i, name, length)
end

"""
Returns the value of the species-save property. By default, concentrations of aqueous species are not saved. Setting the species-save property to true allows aqueous species concentrations to be retrieved with RM_GetSpeciesConcentrations, and solution compositions to be set with RM_SpeciesConcentrations2Module.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0, species are not saved; 1, species are saved; negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_GetSpeciesName, RM_SpeciesConcentrations2Module, RM_SetSpeciesSaveOn.
C Example:
 save_on = RM_GetSpeciesSaveOn(id);
if (save_on .ne. 0)
{
  fprintf(stderr, "Reaction module is saving species concentrations\n");
}
else
{
  fprintf(stderr, "Reaction module is not saving species concentrations\n");
}
"""
function RM_GetSpeciesSaveOn(id::Int)
  IRM_RESULT=ccall((:RM_GetSpeciesSaveOn, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Transfers the charge of each aqueous species to the array argument (z). This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true.

Parameters
id	The instance id returned from RM_Create.
z	Array that receives the charge for each aqueous species. Dimension of the array is nspecies, where nspecies is is the number of aqueous species (RM_GetSpeciesCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesName, RM_SpeciesConcentrations2Module, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
z = (double*) malloc((size_t) (nspecies * sizeof(double)));
status = RM_GetSpeciesZ(id, z);
"""
function RM_GetSpeciesZ(id::Int, z::Array{Float64,1})
  IRM_RESULT=ccall((:RM_GetSpeciesZ, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, z)
end

"""
Returns an array with the starting cell numbers from the range of cell numbers assigned to each worker.

Parameters
id	The instance id returned from RM_Create.
sc	Array to receive the starting cell numbers. Dimension of the array is the number of threads (OpenMP) or the number of processes (MPI).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_Create, RM_GetEndCell, RM_GetMpiTasks, RM_GetThreadCount.
C Example:
 n = RM_GetThreadCount(id) * RM_GetMpiTasks(id);
sc = (int*) malloc((size_t) (n * sizeof(int)));
status = RM_GetStartCell(id, sc);
"""
function RM_GetStartCell(id::Int, sc::Array{Int32,1})
  IRM_RESULT=ccall((:RM_GetStartCell, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cint}), id, sc)
end

"""
Returns the number of threads, which is equal to the number of workers used to run in parallel with OPENMP. For the OPENMP version, the number of threads is set implicitly or explicitly with RM_Create. For the MPI version, the number of threads is always one for each process.

Parameters
id	The instance id returned from RM_Create.
Return values
The	number of threads, negative is failure (See RM_DecodeError).
See also
RM_GetMpiTasks.
C Example:
 sprintf(str1, "Number of threads: %d\n", RM_GetThreadCount(id));
status = RM_OutputMessage(id, str1);
"""
function RM_GetThreadCount(id::Int)
  n_threads=ccall((:RM_GetThreadCount, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the current simulation time in seconds. The reaction module does not change the time value, so the returned value is equal to the default (0.0) or the last time set by RM_SetTime.

Parameters
id	The instance id returned from RM_Create.
Return values
The	current simulation time in seconds.
See also
RM_GetTimeConversion, RM_GetTimeStep, RM_SetTime, RM_SetTimeConversion, RM_SetTimeStep.
C Example:
 sprintf(str, "%s%10.1f%s", "Beginning reaction calculation ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
"""
function RM_GetTime(id::Int)
  t_simulation=ccall((:RM_GetTime, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns a multiplier to convert time from seconds to another unit, as specified by the user. The reaction module uses seconds as the time unit. The user can set a conversion factor (RM_SetTimeConversion) and retrieve it with RM_GetTimeConversion. The reaction module only uses the conversion factor when printing the long version of cell chemistry (RM_SetPrintChemistryOn), which is rare. Default conversion factor is 1.0.

Parameters
id	The instance id returned from RM_Create.
Return values
Multiplier	to convert seconds to another time unit.
See also
RM_GetTime, RM_GetTimeStep, RM_SetTime, RM_SetTimeConversion, RM_SetTimeStep.
C Example:
 sprintf(str, "%s%10.1f%s", "Beginning reaction calculation ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
"""
function RM_GetTimeConversion(id::Int)
  Multiplier=ccall((:RM_GetTimeConversion, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Returns the current simulation time step in seconds. This is the time over which kinetic reactions are integrated in a call to RM_RunCells. The reaction module does not change the time step value, so the returned value is equal to the default (0.0) or the last time step set by RM_SetTimeStep.

Parameters
id	The instance id returned from RM_Create.
Return values
The	current simulation time step in seconds.
See also
RM_GetTime, RM_GetTimeConversion, RM_SetTime, RM_SetTimeConversion, RM_SetTimeStep.
C Example:
 sprintf(str, "%s%10.1f%s", "          Time step                  ",
        RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
"""
function RM_GetTimeStep(id::Int)
  t_step=ccall((:RM_GetTimeStep, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Fills an array (c) with concentrations from solutions in the InitialPhreeqc instance. The method is used to obtain concentrations for boundary conditions. If a negative value is used for a cell in boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance will be used for that cell. Concentrations may be a mixture of two solutions, boundary_solution1 and boundary_solution2, with a mixing fraction for boundary_solution1 1 of fraction1 and mixing fraction for boundary_solution2 of (1 - fraction1). A negative value for boundary_solution2 implies no mixing, and the associated value for fraction1 is ignored. If boundary_solution2 and fraction1 are NULL, no mixing is used; concentrations are derived from boundary_solution1 only.

Parameters
id	The instance id returned from RM_Create.
c	Array of concentrations extracted from the InitialPhreeqc instance. The dimension of c is equivalent to Fortran allocation (n_boundary, ncomp), where ncomp is the number of components returned from RM_FindComponents or RM_GetComponentCount.
n_boundary	The number of boundary condition solutions that need to be filled.
boundary_solution1	Array of solution index numbers that refer to solutions in the InitialPhreeqc instance. Size is n_boundary.
boundary_solution2	Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance and are defined to mix with boundary_solution1. Size is n_boundary. May be NULL in C.
fraction1	Fraction of boundary_solution1 that mixes with (1-fraction1) of boundary_solution2. Size is (n_boundary). May be NULL in C.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetComponentCount.
C Example:
 nbound = 1;
bc1 = (int*) malloc((size_t) (nbound * sizeof(int)));
bc2 = (int*) malloc((size_t) (nbound * sizeof(int)));
bc_f1 = (double*) malloc((size_t) (nbound * sizeof(double)));
bc_conc = (double*) malloc((size_t) (ncomps * nbound * sizeof(double)));
for (i = 0; i < nbound; i++)
{
  bc1[i]          = 0;       // Solution 0 from InitialPhreeqc instance
  bc2[i]          = -1;      // no bc2 solution for mixing
  bc_f1[i]        = 1.0;     // mixing fraction for bc1
}
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);
"""
function RM_InitialPhreeqc2Concentrations(id::Int, c::Array{Float64,1},
  n_boundary::Int, boundary_solution1::Array{Int,1},
  boundary_solution2::Array{Int,1}, fraction1::Array{Float64,1})
  boundary_solution1=convert(Array{Int32}, boundary_solution1)
  boundary_solution2=convert(Array{Int32}, boundary_solution2)
  IRM_RESULT=ccall((:RM_InitialPhreeqc2Concentrations, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
  id,c,n_boundary,boundary_solution1,boundary_solution2,fraction1)
end

"""
Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers, possibly with mixing. In its simplest form, initial_conditions1 is used to select initial conditions, including solutions and reactants, for each cell of the model, without mixing. Initial_conditions1 is dimensioned (nxyz, 7), where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order: (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS. In C, initial_solution1[3*nxyz + 99] = 2, indicates that cell 99 (0 based) contains the SURFACE definition with user number 2 that has been defined in the InitialPhreeqc instance (either by RM_RunFile or RM_RunString).

It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing, initials_conditions2 contains numbers for a second entity that mixes with the entity defined in initial_conditions1. Fraction1 contains the mixing fraction for initial_conditions1, whereas (1 - fraction1) is the mixing fraction for initial_conditions2. In C, initial_solution1[3*nxyz + 99] = 2, initial_solution2[3*nxyz + 99] = 3, fraction1[3*nxyz + 99] = 0.25 indicates that cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface compositions have been defined in the InitialPhreeqc instance. If the user number in initial_conditions2 is negative, no mixing occurs. If initials_conditions2 and fraction1 are NULL, no mixing is used, and initial conditions are derived solely from initials_conditions1.

Parameters
id	The instance id returned from RM_Create.
initial_conditions1	Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance. Size is (nxyz,7). The order of definitions is given above. Negative values are ignored, resulting in no definition of that entity for that cell.
initial_conditions2	Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance. Nonnegative values of initial_conditions2 result in mixing with the entities defined in initial_conditions1. Negative values result in no mixing. Size is (nxyz,7). The order of definitions is given above. May be NULL in C; setting to NULL results in no mixing.
fraction1	Fraction of initial_conditions1 that mixes with (1-fraction1) of initial_conditions2. Size is (nxyz,7). The order of definitions is given above. May be NULL in C; setting to NULL results in no mixing.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqcCell2Module.
C Example:
 ic1 = (int*) malloc((size_t) (7 * nxyz * sizeof(int)));
ic2 = (int*) malloc((size_t) (7 * nxyz * sizeof(int)));
f1 = (double*) malloc((size_t) (7 * nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
  ic1[i]          = 1;       // Solution 1
  ic1[nxyz + i]   = -1;      // Equilibrium phases none
  ic1[2*nxyz + i] = 1;       // Exchange 1
  ic1[3*nxyz + i] = -1;      // Surface none
  ic1[4*nxyz + i] = -1;      // Gas phase none
  ic1[5*nxyz + i] = -1;      // Solid solutions none
  ic1[6*nxyz + i] = -1;      // Kinetics none

  ic2[i]          = -1;      // Solution none
  ic2[nxyz + i]   = -1;      // Equilibrium phases none
  ic2[2*nxyz + i] = -1;      // Exchange none
  ic2[3*nxyz + i] = -1;      // Surface none
  ic2[4*nxyz + i] = -1;      // Gas phase none
  ic2[5*nxyz + i] = -1;      // Solid solutions none
  ic2[6*nxyz + i] = -1;      // Kinetics none

  f1[i]          = 1.0;      // Mixing fraction ic1 Solution
  f1[nxyz + i]   = 1.0;      // Mixing fraction ic1 Equilibrium phases
  f1[2*nxyz + i] = 1.0;      // Mixing fraction ic1 Exchange 1
  f1[3*nxyz + i] = 1.0;      // Mixing fraction ic1 Surface
  f1[4*nxyz + i] = 1.0;      // Mixing fraction ic1 Gas phase
  f1[5*nxyz + i] = 1.0;      // Mixing fraction ic1 Solid solutions
  f1[6*nxyz + i] = 1.0;      // Mixing fraction ic1 Kinetics
}
status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1);
// No mixing is defined, so the following is equivalent
status = RM_InitialPhreeqc2Module(id, ic1, NULL, NULL);
"""
function RM_InitialPhreeqc2Module(id::Int, initial_conditions1::Array{Int,1},
  initial_conditions2::Array{Int,1}, fraction1::Array{Float64,1})
  initial_conditions1=convert(Array{Int32}, initial_conditions1)
  initial_conditions2=convert(Array{Int32}, initial_conditions2)
  IRM_RESULT=ccall((:RM_InitialPhreeqc2Module, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}), id,initial_conditions1,
  initial_conditions2,fraction1)
end

"""
Fills an array (species_c) with aqueous species concentrations from solutions in the InitialPhreeqc instance. This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value is used for a cell in boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance will be used for that cell. Concentrations may be a mixture of two solutions, boundary_solution1 and boundary_solution2, with a mixing fraction for boundary_solution1 1 of fraction1 and mixing fraction for boundary_solution2 of (1 - fraction1). A negative value for boundary_solution2 implies no mixing, and the associated value for fraction1 is ignored. If boundary_solution2 and fraction1 are NULL, no mixing is used; concentrations are derived from boundary_solution1 only.

Parameters
id	The instance id returned from RM_Create.
species_c	Array of aqueous concentrations extracted from the InitialPhreeqc instance. The dimension of species_c is equivalent to Fortran allocation (n_boundary, nspecies), where nspecies is the number of aqueous species returned from RM_GetSpeciesCount.
n_boundary	The number of boundary condition solutions that need to be filled.
boundary_solution1	Array of solution index numbers that refer to solutions in the InitialPhreeqc instance. Size is n_boundary.
boundary_solution2	Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance and are defined to mix with boundary_solution1. Size is n_boundary. May be NULL in C.
fraction1	Fraction of boundary_solution1 that mixes with (1-fraction1) of boundary_solution2. Size is n_boundary. May be NULL in C.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesCount, RM_SetSpeciesSaveOn.
C Example:
 nbound = 1;
nspecies = RM_GetSpeciesCount(id);
bc1 = (int*) malloc((size_t) (nbound * sizeof(int)));
bc2 = (int*) malloc((size_t) (nbound * sizeof(int)));
bc_f1 = (double*) malloc((size_t) (nbound * sizeof(double)));
bc_conc = (double*) malloc((size_t) (nspecies * nbound * sizeof(double)));
for (i = 0; i < nbound; i++)
{
  bc1[i]          = 0;       // Solution 0 from InitialPhreeqc instance
  bc2[i]          = -1;      // no bc2 solution for mixing
  bc_f1[i]        = 1.0;     // mixing fraction for bc1
}
status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);
"""
function RM_InitialPhreeqc2SpeciesConcentrations(id::Int,
  species_c::Array{Float64,1}, n_boundary::Int, boundary_solution1::Array{Int,1},
  boundary_solution2::Array{Int,1}, fraction1::Array{Float64,1})
  boundary_solution1=convert(Array{Int32}, boundary_solution1)
  boundary_solution2=convert(Array{Int32}, boundary_solution2)
  IRM_RESULT=ccall((:RM_InitialPhreeqc2SpeciesConcentrations, Lib_PhreeqcRM_path),
  Cint,  (Cint,Ptr{Cdouble},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
  id, species_c, n_boundary, boundary_solution1, boundary_solution2, fraction1)
end

"""
A cell numbered n in the InitialPhreeqc instance is selected to populate a series of cells. All reactants with the number n are transferred along with the solution. If MIX n exists, it is used for the definition of the solution. If n is negative, n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance. All reactants for each cell in the list module_numbers are removed before the cell definition is copied from the InitialPhreeqc instance to the workers.

Parameters
id	The instance id returned from RM_Create.
n	Cell number refers to a solution or MIX and associated reactants in the InitialPhreeqc instance. A negative number indicates the largest solution or MIX number in the InitialPhreeqc instance will be used.
module_numbers	A list of cell numbers in the user's grid-cell numbering system that will be populated with cell n from the InitialPhreeqc instance.
dim_module_numbers	The number of cell numbers in the module_numbers list.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module.
C Example:
 module_cells = (int*) malloc((size_t) (2 * sizeof(int)));
module_cells[0] = 18;
module_cells[1] = 19;
// n will be the largest SOLUTION number in InitialPhreeqc instance
// copies solution and reactants to cells 18 and 19
status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2);
"""
function RM_InitialPhreeqcCell2Module(id::Int, n::Int,
  module_numbers::Array{Int,1}, dim_module_numbers::Int)
  module_numbers=convert(Array{Int32}, module_numbers)
  IRM_RESULT=ccall((:RM_InitialPhreeqcCell2Module, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Ptr{Cint},Cint), id, n, module_numbers, dim_module_numbers)
end

"""
Load a database for all IPhreeqc instances–workers, InitialPhreeqc, and Utility. All definitions of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.

Parameters
id	The instance id returned from RM_Create.
db_name	String containing the database name.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_Create.
C Example:
 status = RM_LoadDatabase(id, "phreeqc.dat");
"""
function RM_LoadDatabase(id::Int, db_name::AbstractString)
  IRM_RESULT=ccall((:RM_LoadDatabase, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, db_name)
end

"""
Print a message to the log file.

Parameters
id	The instance id returned from RM_Create.
str	String to be printed.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_OpenFiles, RM_ErrorMessage, RM_OutputMessage, RM_ScreenMessage, RM_WarningMessage.
C Example:
 sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
"""
function RM_LogMessage(id::Int, str::AbstractString)
  IRM_RESULT=ccall((:RM_LogMessage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, str)
end

"""
MPI only. Workers (processes with RM_GetMpiMyself > 0) must call RM_MpiWorker to be able to respond to messages from the root to accept data, perform calculations, and (or) return data. RM_MpiWorker contains a loop that reads a message from root, performs a task, and waits for another message from root. RM_SetConcentrations, RM_RunCells, and RM_GetConcentrations are examples of methods that send a message from root to get the workers to perform a task. The workers will respond to all methods that are designated "workers must be in the loop of RM_MpiWorker" in the MPI section of the method documentation. The workers will continue to respond to messages from root until root calls RM_MpiWorkerBreak.

(Advanced) The list of tasks that the workers perform can be extended by using RM_SetMpiWorkerCallback. It is then possible to use the MPI processes to perform other developer-defined tasks, such as transport calculations, without exiting from the RM_MpiWorker loop. Alternatively, root calls RM_MpiWorkerBreak to allow the workers to continue past a call to RM_MpiWorker. The workers perform developer-defined calculations, and then RM_MpiWorker is called again to respond to requests from root to perform reaction-module tasks.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError). RM_MpiWorker returns a value only when RM_MpiWorkerBreak is called by root.
See also
RM_MpiWorkerBreak, RM_SetMpiWorkerCallback, RM_SetMpiWorkerCallbackCookie.
C Example:
 status = RM_MpiWorker(id);
"""
function RM_MpiWorker(id::Int)
  IRM_RESULT=ccall((:RM_MpiWorker, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
MPI only. This method is called by root to force workers (processes with RM_GetMpiMyself > 0) to return from a call to RM_MpiWorker. RM_MpiWorker contains a loop that reads a message from root, performs a task, and waits for another message from root. The workers respond to all methods that are designated "workers must be in the loop of RM_MpiWorker" in the MPI section of the method documentation. The workers will continue to respond to messages from root until root calls RM_MpiWorkerBreak.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_MpiWorker, RM_SetMpiWorkerCallback, RM_SetMpiWorkerCallbackCookie.
C Example:
 status = RM_MpiWorkerBreak(id);
"""
function RM_MpiWorkerBreak(id::Int)
  IRM_RESULT=ccall((:RM_MpiWorkerBreak, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt based on the prefix defined by RM_SetFilePrefix.

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetFilePrefix, RM_GetFilePrefix, RM_CloseFiles, RM_ErrorMessage, RM_LogMessage, RM_OutputMessage, RM_WarningMessage.
C Example:
 status = RM_SetFilePrefix(id, "Advect_c");
status = RM_OpenFiles(id);
"""
function RM_OpenFiles(id::Int)
  IRM_RESULT=ccall((:RM_OpenFiles, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Print a message to the output file.

Parameters
id	The instance id returned from RM_Create.
str	String to be printed.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_ErrorMessage, RM_LogMessage, RM_ScreenMessage, RM_WarningMessage.
C Example:
 sprintf(str1, "Number of threads:                                %d\n", RM_GetThreadCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of MPI processes:                          %d\n", RM_GetMpiTasks(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "MPI task number:                                  %d\n", RM_GetMpiMyself(id));
status = RM_OutputMessage(id, str1);
status = RM_GetFilePrefix(id, str, 100);
sprintf(str1, "File prefix:                                      %s\n", str);
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of grid cells in the user's model:         %d\n", RM_GetGridCellCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of chemistry cells in the reaction module: %d\n", RM_GetChemistryCellCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of components for transport:               %d\n", RM_GetComponentCount(id));
status = RM_OutputMessage(id, str1);
"""
function RM_OutputMessage(id::Int, str::AbstractString)
  IRM_RESULT=ccall((:RM_OutputMessage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, str)
end

"""
Runs a reaction step for all of the cells in the reaction module. Normally, tranport concentrations are transferred to the reaction cells (RM_SetConcentrations) before reaction calculations are run. The length of time over which kinetic reactions are integrated is set by RM_SetTimeStep. Other properties that may need to be updated as a result of the transport calculations include porosity (RM_SetPorosity), saturation (RM_SetSaturation), temperature (RM_SetTemperature), and pressure (RM_SetPressure).

Parameters
id	The instance id returned from RM_Create.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetConcentrations, RM_SetPorosity, RM_SetPressure, RM_SetSaturation, RM_SetTemperature, RM_SetTimeStep.
C Example:
 status = RM_SetPorosity(id, por);              // If porosity changes
status = RM_SetSaturation(id, sat);            // If saturation changes
status = RM_SetTemperature(id, temperature);   // If temperature changes
status = RM_SetPressure(id, pressure);         // If pressure changes
status = RM_SetConcentrations(id, c);          // Transported concentrations
status = RM_SetTimeStep(id, time_step);        // Time step for kinetic reactions
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);          // Concentrations after reaction
status = RM_GetDensity(id, density);           // Density after reaction
status = RM_GetSolutionVolume(id, volume);     // Solution volume after reaction
"""
function RM_RunCells(id::Int)
  IRM_RESULT=ccall((:RM_RunCells, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""
Run a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run the file–the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input files that modify the thermodynamic database should be run by all three sets of instances. Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to be run by the workers. Files that contain initial conditions or boundary conditions should be run by the InitialPhreeqc instance.

Parameters
id	The instance id returned from RM_Create.
workers	1, the workers will run the file; 0, the workers will not run the file.
initial_phreeqc	1, the InitialPhreeqc instance will run the file; 0, the InitialPhreeqc will not run the file.
utility	1, the Utility instance will run the file; 0, the Utility instance will not run the file.
chem_name	Name of the file to run.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_RunString.
C Example:
 status = RM_RunFile(id, 1, 1, 1, "advect.pqi");
"""
function RM_RunFile(id::Int, workers::Int, initial_phreeqc::Int,
  utility::Int, chem_name::AbstractString)
  IRM_RESULT=ccall((:RM_RunFile, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Cint,Cint,Cstring), id, workers, initial_phreeqc, utility, chem_name)
end

"""
Run a PHREEQC input string. The first three arguments determine which IPhreeqc instances will run the string–the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input strings that modify the thermodynamic database should be run by all three sets of instances. Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to be run by the workers. Strings that contain initial conditions or boundary conditions should be run by the InitialPhreeqc instance.

Parameters
id	The instance id returned from RM_Create.
workers	1, the workers will run the string; 0, the workers will not run the string.
initial_phreeqc	1, the InitialPhreeqc instance will run the string; 0, the InitialPhreeqc will not run the string.
utility	1, the Utility instance will run the string; 0, the Utility instance will not run the string.
input_string	String containing PHREEQC input.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_RunFile.
C Example:
 strcpy(str, "DELETE; -all");
status = RM_RunString(id, 1, 0, 1, str);    // workers, initial_phreeqc, utility
"""
function RM_RunString(id::Int, workers::Int, initial_phreeqc::Int,
  utility::Int, input_string::AbstractString)
  IRM_RESULT=ccall((:RM_RunString, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Cint,Cint,Cstring), id, workers, initial_phreeqc,
  utility, input_string)
end

"""
Print message to the screen.

Parameters
id	The instance id returned from RM_Create.
str	String to be printed.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_ErrorMessage, RM_LogMessage, RM_OutputMessage, RM_WarningMessage.
C Example:
 sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ",
        time * RM_GetTimeConversion(id), " days\n");
status = RM_ScreenMessage(id, str);
"""
function RM_ScreenMessage(id::Int, str::AbstractString)
  IRM_RESULT=ccall((:RM_ScreenMessage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, str)
end

"""
Select whether to include H2O in the component list. The concentrations of H and O must be known accurately (8 to 10 significant digits) for the numerical method of PHREEQC to produce accurate pH and pe values. Because most of the H and O are in the water species, it may be more robust (require less accuracy in transport) to transport the excess H and O (the H and O not in water) and water. The default setting (true) is to include water, excess H, and excess O as components. A setting of false will include total H and total O as components. RM_SetComponentH2O must be called before RM_FindComponents.

Parameters
id	The instance id returned from RM_Create.
tf	0, total H and O are included in the component list; 1, excess H, excess O, and water are included in the component list.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents.
C Example:
 status = RM_SetComponentH2O(id, 0);
"""
function RM_SetComponentH2O(id::Int, tf::Int)
  IRM_RESULT=ccall((:RM_SetComponentH2O, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, tf)
end

"""
Use the vector of concentrations (c) to set the moles of components in each reaction cell. The volume of water in a cell is the product of porosity (RM_SetPorosity), saturation (RM_SetSaturation), and reference volume (RM_SetRepresentativeVolume). The moles of each component are determined by the volume of water and per liter concentrations. If concentration units (RM_SetUnitsSolution) are mass fraction, the density (as specified by RM_SetDensity) is used to convert from mass fraction to per mass per liter.

Parameters
id	The instance id returned from RM_Create.
c	Array of component concentrations. Size of array is equivalent to Fortran (nxyz, ncomps), where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount), and ncomps is the number of components as determined by RM_FindComponents or RM_GetComponentCount.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetDensity, RM_SetPorosity, RM_SetRepresentativeVolume, RM_SetSaturation, RM_SetUnitsSolution.
C Example:
 c = (double*) malloc((size_t) (ncomps * nxyz * sizeof(double)));
...
advect_c(c, bc_conc, ncomps, nxyz, nbound);
status = RM_SetPoosity(id, por);               // If porosity changes
status = RM_SetSaturation(id, sat);            // If saturation changes
status = RM_SetTemperature(id, temperature);   // If temperature changes
status = RM_SetPressure(id, pressure);         // If pressure changes
status = RM_SetConcentrations(id, c);          // Transported concentrations
status = RM_SetTimeStep(id, time_step);        // Time step for kinetic reactions
status = RM_SetTime(id, time);                 // Current time
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);          // Concentrations after reaction
status = RM_GetDensity(id, density);           // Density after reaction
status = RM_GetSolutionVolume(id, volume);     // Solution volume after reaction
"""
function RM_SetConcentrations(id::Int, c::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetConcentrations, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, c)
end

"""
Select the current selected output by user number. The user may define multiple SELECTED_OUTPUT data blocks for the workers. A user number is specified for each data block. The value of the argument n_user selects which of the SELECTED_OUTPUT definitions will be used for selected-output operations.

Parameters
id	The instance id returned from RM_Create.
n_user	User number of the SELECTED_OUTPUT data block that is to be used.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetNthSelectedOutputUserNumber, RM_GetSelectedOutput, RM_GetSelectedOutputColumnCount, RM_GetSelectedOutputCount, RM_GetSelectedOutputRowCount, RM_GetSelectedOutputHeading, RM_SetSelectedOutputOn.
C Example:
 for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double*) malloc((size_t) (col* nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
"""
function RM_SetCurrentSelectedOutputUserNumber(id::Int, n_user::Int)
  IRM_RESULT=ccall((:RM_SetCurrentSelectedOutputUserNumber, Lib_PhreeqcRM_path),
  Cint,  (Cint,Cint), id, n_user)
end

"""
Set the density for each reaction cell. These density values are used when converting from transported mass fraction concentrations (RM_SetUnitsSolution) to produce per liter concentrations during a call to RM_SetConcentrations. They are also used when converting from module concentrations to transport concentrations of mass fraction (RM_GetConcentrations), if RM_UseSolutionDensityVolume is set to false.

Parameters
id	The instance id returned from RM_Create.
density	Array of densities. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetConcentrations, RM_SetConcentrations, RM_SetUnitsSolution, RM_UseSolutionDensityVolume.
C Example:
 density = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
    density[i] = 1.0;
}
status = RM_SetDensity(id, density);
"""
function RM_SetDensity(id::Int, density::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetDensity, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, density)
end

"""
Set the name of the dump file. It is the name used by RM_DumpModule.

Parameters
id	The instance id returned from RM_Create.
dump_name	Name of dump file.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_DumpModule.
C Example:
 status = RM_SetDumpFileName(id, "advection_c.dmp");
dump_on = 1;
append = 0;
status = RM_DumpModule(id, dump_on, append);
"""
function RM_SetDumpFileName(id::Int, dump_name::AbstractString)
  IRM_RESULT=ccall((:RM_SetDumpFileName, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, dump_name)
end

"""
Set the action to be taken when the reaction module encounters an error. Options are 0, return to calling program with an error return code (default); 1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or 2, attempt to exit gracefully.

Parameters
id	The instance id returned from RM_Create.
mode	Error handling mode: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 id = RM_Create(nxyz, nthreads);
status = RM_SetErrorHandlerMode(id, 2);
"""
function RM_SetErrorHandlerMode(id::Int, mode::Int)
  IRM_RESULT=ccall((:RM_SetErrorHandlerMode, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, mode)
end

"""
Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files. These files are opened by RM_OpenFiles.

Parameters
id	The instance id returned from RM_Create.
prefix	Prefix used when opening the output and log files.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_OpenFiles, RM_CloseFiles.
C Example:
 status = RM_SetFilePrefix(id, "Advect_c");
status = RM_OpenFiles(id);
"""
function RM_SetFilePrefix(id::Int, prefix::AbstractString)
  IRM_RESULT=ccall((:RM_SetFilePrefix, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, prefix)
end

# MPI functions are not wrapped!
"""
Sets the property for partitioning solids between the saturated and unsaturated parts of a partially saturated cell.

The option is intended to be used by saturated-only flow codes that allow a variable water table. The value has meaning only when saturations less than 1.0 are encountered. The partially saturated cells may have a small water-to-rock ratio that causes reactions to proceed differently relative to fully saturated cells. By setting RM_SetPartitionUZSolids to true, the amounts of solids and gases are partioned according to the saturation. If a cell has a saturation of 0.5, then the water interacts with only half of the solids and gases; the other half is unreactive until the water table rises. As the saturation in a cell varies, solids and gases are transferred between the saturated and unsaturated (unreactive) reservoirs of the cell. Unsaturated-zone flow and transport codes will probably use the default (false), which assumes all gases and solids are reactive regardless of saturation.

Parameters
id	The instance id returned from RM_Create.
tf	True, the fraction of solids and gases available for reaction is equal to the saturation; False (default), all solids and gases are reactive regardless of saturation.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
C Example:
 status = RM_SetPartitionUZSolids(id, 0);
"""
function RM_SetPartitionUZSolids(id::Int, tf::Int)
  IRM_RESULT=ccall((:RM_SetPartitionUZSolids, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, tf)
end

"""
Set the porosity for each reaction cell. The volume of water in a reaction cell is the product of the porosity, the saturation (RM_SetSaturation), and the representative volume (RM_SetRepresentativeVolume).

Parameters
id	The instance id returned from RM_Create.
por	Array of porosities, unitless. Default is 0.1. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetSaturation, RM_SetRepresentativeVolume, RM_SetSaturation.
C Example:
 por = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) por[i] = 0.2;
status = RM_SetPorosity(id, por);
"""
function RM_SetPorosity(id::Int, por::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetPorosity, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, por)
end

"""
Set the pressure for each reaction cell. Pressure effects are considered only in three of the databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.

Parameters
id	The instance id returned from RM_Create.
p	Array of pressures, in atm. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetTemperature.
C Example:
 pressure = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) pressure[i] = 2.0;
status = RM_SetPressure(id, pressure);
"""
function RM_SetPressure(id::Int, p::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetPressure, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, p)
end

"""
Enable or disable detailed output for each reaction cell. Printing for a cell will occur only when the printing is enabled with RM_SetPrintChemistryOn and the cell_mask value is 1.

Parameters
id	The instance id returned from RM_Create.
cell_mask	Array of integers. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount). A value of 0 will disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetPrintChemistryOn.
C Example:
 print_chemistry_mask = (int*) malloc((size_t) (nxyz * sizeof(int)));
for (i = 0; i < nxyz/2; i++)
{
  print_chemistry_mask[i] = 1;
  print_chemistry_mask[i + nxyz/2] = 0;
}
status = RM_SetPrintChemistryMask(id, print_chemistry_mask);
"""
function RM_SetPrintChemistryMask(id::Int, cell_mask::Array{Int,1})
  cell_mask=convert(Array{Int32}, cell_mask)
  IRM_RESULT=ccall((:RM_SetPrintChemistryMask, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cint}), id, cell_mask)
end

"""
Setting to enable or disable printing detailed output from reaction calculations to the output file for a set of cells defined by RM_SetPrintChemistryMask. The detailed output prints all of the output typical of a PHREEQC reaction calculation, which includes solution descriptions and the compositions of all other reactants. The output can be several hundred lines per cell, which can lead to a very large output file (prefix.chem.txt, RM_OpenFiles). For the worker instances, the output can be limited to a set of cells (RM_SetPrintChemistryMask) and, in general, the amount of information printed can be limited by use of options in the PRINT data block of PHREEQC (applied by using RM_RunFile or RM_RunString). Printing the detailed output for the workers is generally used only for debugging, and PhreeqcRM will run significantly faster when printing detailed output for the workers is disabled.

Parameters
id	The instance id returned from RM_Create.
workers	0, disable detailed printing in the worker instances, 1, enable detailed printing in the worker instances.
initial_phreeqc	0, disable detailed printing in the InitialPhreeqc instance, 1, enable detailed printing in the InitialPhreeqc instances.
utility	0, disable detailed printing in the Utility instance, 1, enable detailed printing in the Utility instance.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetPrintChemistryMask.
C Example:
 status = RM_SetPrintChemistryOn(id, 0, 1, 0); // workers, initial_phreeqc, utility
"""
function RM_SetPrintChemistryOn(id::Int, workers::Int,
  initial_phreeqc::Int, utility::Int)
  IRM_RESULT=ccall((:RM_SetPrintChemistryOn, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint,Cint,Cint), id, workers, initial_phreeqc, utility)
end

"""
Set the load-balancing algorithm. PhreeqcRM attempts to rebalance the load of each thread or process such that each thread or process takes the same amount of time to run its part of a RM_RunCells calculation. Two algorithms are available; one uses individual times for each cell and accounts for cells that were not run because saturation was zero (default), and the other assigns an average time to all cells. The methods are similar, but limited testing indicates the default method performs better.

Parameters
id	The instance id returned from RM_Create.
method	0, indicates average times are used in rebalancing; 1 indicates individual cell times are used in rebalancing (default).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetRebalanceFraction.
C Example:
 status = RM_SetRebalanceByCell(id, 1);
"""
function RM_SetRebalanceByCell(id::Int, method::Int)
  IRM_RESULT=ccall((:RM_SetRebalanceByCell, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, method)
end

"""
Sets the fraction of cells that are transferred among threads or processes when rebalancing. PhreeqcRM attempts to rebalance the load of each thread or process such that each thread or process takes the same amount of time to run its part of a RM_RunCells calculation. The rebalancing transfers cell calculations among threads or processes to try to achieve an optimum balance. RM_SetRebalanceFraction adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to determine the actual number of cell transfers. A value of zero eliminates load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell distribution and avoid possible oscillations when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration. Default is 0.5.

Parameters
id	The instance id returned from RM_Create.
f	Fraction from 0.0 to 1.0.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetRebalanceByCell.
C Example:
 status = RM_SetRebalanceFraction(id, 0.5);
"""
function RM_SetRebalanceFraction(id::Int, f::Float64)
  IRM_RESULT=ccall((:RM_SetRebalanceFraction, Lib_PhreeqcRM_path), Cint,
  (Cint,Cdouble), id, f)
end

"""
Set the representative volume of each reaction cell. By default the representative volume of each reaction cell is 1 liter. The volume of water in a reaction cell is determined by the procuct of the representative volume, the porosity (RM_SetPorosity), and the saturation (RM_SetSaturation). The numerical method of PHREEQC is more robust if the water volume for a reaction cell is within a couple orders of magnitude of 1.0. Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes) may cause non-convergence of the numerical method. In these cases, a larger representative volume may help. Note that increasing the representative volume also increases the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers, and others), which are defined as moles per representative volume.

Parameters
id	The instance id returned from RM_Create.
rv	Vector of representative volumes, in liters. Default is 1.0 liter. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetPorosity, RM_SetSaturation.
C Example:
 double * rv;
rv = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) rv[i] = 1.0;
status = RM_SetRepresentativeVolume(id, rv);
"""
function RM_SetRepresentativeVolume(id::Int, rv::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetRepresentativeVolume, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, rv)
end

"""
Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1. The volume of water in a cell is the product of porosity (RM_SetPorosity), saturation (RM_SetSaturation), and representative volume (RM_SetRepresentativeVolume). As a result of a reaction calculation, solution properties (density and volume) will change; the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume data to calculate these changes. The methods RM_GetDensity, RM_GetSolutionVolume, and RM_GetSaturation can be used to account for these changes in the succeeding transport calculation. RM_SetRepresentativeVolume should be called before initial conditions are defined for the reaction cells.

Parameters
id	The instance id returned from RM_Create.
sat	Array of saturations, unitless. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetDensity, RM_GetSaturation, RM_GetSolutionVolume, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 sat = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) sat[i] = 1.0;
status = RM_SetSaturation(id, sat);
"""
function RM_SetSaturation(id::Int, sat::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetSaturation, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, sat)
end

"""
Set the property that controls whether messages are written to the screen. Messages include information about rebalancing during RM_RunCells, and any messages written with RM_ScreenMessage.

Parameters
id	The instance id returned from RM_Create.
tf	1, enable screen messages; 0, disable screen messages. Default is 1.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_RunCells, RM_ScreenMessage.
C Example:
 status = RM_SetScreenOn(id, 1);
"""
function RM_SetScreenOn(id::Int, tf::Int)
  IRM_RESULT=ccall((:RM_SetScreenOn, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id,tf)
end

"""
Setting determines whether selected-output results are available to be retrieved with RM_GetSelectedOutput. 1 indicates that selected-output results will be accumulated during RM_RunCells and can be retrieved with RM_GetSelectedOutput; 0 indicates that selected-output results will not be accumulated during RM_RunCells.

Parameters
id	The instance id returned from RM_Create.
selected_output	0, disable selected output; 1, enable selected output.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_GetSelectedOutput, RM_SetPrintChemistryOn.
C Example:
 status = RM_SetSelectedOutputOn(id, 1);       // enable selected output
"""
function RM_SetSelectedOutputOn(id::Int, selected_output::Int)
  IRM_RESULT=ccall((:RM_SetSelectedOutputOn, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, selected_output)
end

"""
Sets the value of the species-save property. This method enables use of PhreeqcRM with multicomponent-diffusion transport calculations. By default, concentrations of aqueous species are not saved. Setting the species-save property to 1 allows aqueous species concentrations to be retrieved with RM_GetSpeciesConcentrations, and solution compositions to be set with RM_SpeciesConcentrations2Module. RM_SetSpeciesSaveOn must be called before calls to RM_FindComponents.

Parameters
id	The instance id returned from RM_Create.
save_on	0, indicates species concentrations are not saved; 1, indicates species concentrations are saved.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesSaveOn, RM_GetSpeciesZ, RM_GetSpeciesName, RM_SpeciesConcentrations2Module.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
"""
function RM_SetSpeciesSaveOn(id::Int, save_on::Int)
  IRM_RESULT=ccall((:RM_SetSpeciesSaveOn, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, save_on)
end

"""
Set the temperature for each reaction cell. If RM_SetTemperature is not called, worker solutions will have temperatures as defined by initial conditions (RM_InitialPhreeqc2Module and RM_InitialPhreeqcCell2Module).

Parameters
id	The instance id returned from RM_Create.
t	Array of temperatures, in degrees C. Size of array is nxyz, where nxyz is the number of grid cells in the user's model (RM_GetGridCellCount).
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPressure.
C Example:
 temperature = (double*) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
  temperature[i] = 20.0;
}
status = RM_SetTemperature(id, temperature);
"""
function RM_SetTemperature(id::Int, t::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SetTemperature, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, t)
end

"""
Set current simulation time for the reaction module.

Parameters
id	The instance id returned from RM_Create.
time	Current simulation time, in seconds.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetTimeStep, RM_SetTimeConversion.
C Example:
 status = RM_SetTime(id, time);
"""
function RM_SetTime(id::Int, t::Float64)
  IRM_RESULT=ccall((:RM_SetTime, Lib_PhreeqcRM_path), Cint,
  (Cint,Cdouble), id, t)
end

"""
Set a factor to convert to user time units. Factor times seconds produces user time units.

Parameters
id	The instance id returned from RM_Create.
conv_factor	Factor to convert seconds to user time units.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetTime, RM_SetTimeStep.
C Example:
 status = RM_SetTimeConversion(id, 1.0 / 86400.0); // days
"""
function RM_SetTimeConversion(id::Int, conv_factor::Float64)
  IRM_RESULT=ccall((:RM_SetTimeConversion, Lib_PhreeqcRM_path), Cint,
  (Cint,Cdouble), id, conv_factor)
end

"""
Set current time step for the reaction module. This is the length of time over which kinetic reactions are integrated.

Parameters
id	The instance id returned from RM_Create.
time_step	Current time step, in seconds.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetTime, RM_SetTimeConversion.
C Example:
 time_step = 86400.0;
status = RM_SetTimeStep(id, time_step);
"""
function RM_SetTimeStep(id::Int, time_step::Float64)
  IRM_RESULT=ccall((:RM_SetTimeStep, Lib_PhreeqcRM_path), Cint,
  (Cint,Cdouble), id, time_step)
end

"""
Sets input units for exchangers. In PHREEQC input, exchangers are defined by moles of exchange sites (Mp). RM_SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single EXCHANGE definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of exchangers will be the same regardless of porosity. For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.

Parameters
id	The instance id returned from RM_Create.
option	Units option for exchangers: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 status = RM_SetUnitsExchange(id, 1);
"""
function RM_SetUnitsExchange(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsExchange, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set input units for gas phases. In PHREEQC input, gas phases are defined by moles of component gases (Mp). RM_SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single GAS_PHASE definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of a gas component will be the same regardless of porosity. For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.

Parameters
id	The instance id returned from RM_Create.
option	Units option for gas phases: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 status = RM_SetUnitsGasPhase(id, 1);
"""
function RM_SetUnitsGasPhase(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsGasPhase, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set input units for kinetic reactants.

In PHREEQC input, kinetics are defined by moles of kinetic reactants (Mp). RM_SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single KINETICS definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of kinetic reactants will be the same regardless of porosity. For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.

Note that the volume of water in a cell in the reaction module is equal to the product of porosity (RM_SetPorosity), the saturation (RM_SetSaturation), and representative volume (RM_SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the RATES definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of water, often by calculating the rate of reaction per liter of water and multiplying by the volume of water (Basic function SOLN_VOL).

Rates that depend on surface area of solids, are not dependent on the volume of water. However, it is important to get the correct surface area for the kinetic reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant) can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of reactant (Basic function M) in RATES to obtain the surface area.

Parameters
id	The instance id returned from RM_Create.
option	Units option for kinetic reactants: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume, RM_SetSaturation.
C Example:
 status = RM_SetUnitsKinetics(id, 1);
"""
function RM_SetUnitsKinetics(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsKinetics, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set input units for pure phase assemblages (equilibrium phases). In PHREEQC input, equilibrium phases are defined by moles of each phase (Mp). RM_SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of a mineral will be the same regardless of porosity. For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.

Parameters
id	The instance id returned from RM_Create.
option	Units option for equilibrium phases: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 status = RM_SetUnitsPPassemblage(id, 1);
"""
function RM_SetUnitsPPassemblage(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsPPassemblage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Solution concentration units used by the transport model. Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs. PHREEQC defines solutions by the number of moles of each element in the solution.

To convert from mg/L to moles of element in the representative volume of a reaction cell, mg/L is converted to mol/L and multiplied by the solution volume, which is the product of porosity (RM_SetPorosity), saturation (RM_SetSaturation), and representative volume (RM_SetRepresentativeVolume). To convert from mol/L to moles of element in the representative volume of a reaction cell, mol/L is multiplied by the solution volume. To convert from mass fraction to moles of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density (RM_SetDensity) and multiplied by the solution volume.

To convert from moles of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the solution volume resulting in mol/L, and then converted to mg/L. To convert from moles of element in a cell to mol/L, the number of moles of an element is divided by the solution volume resulting in mol/L. To convert from moles of element in a cell to mass fraction, the number of moles of an element is converted to kg and divided by the total mass of the solution. Two options are available for the volume and mass of solution that are used in converting to transport concentrations: (1) the volume and mass of solution are calculated by PHREEQC, or (2) the volume of solution is the product of porosity (RM_SetPorosity), saturation (RM_SetSaturation), and representative volume (RM_SetRepresentativeVolume), and the mass of solution is volume times density as defined by RM_SetDensity. Which option is used is determined by RM_UseSolutionDensityVolume.

Parameters
id	The instance id returned from RM_Create.
option	Units option for solutions: 1, 2, or 3, default is 1, mg/L.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_SetDensity, RM_SetPorosity, RM_SetRepresentativeVolume, RM_SetSaturation, RM_UseSolutionDensityVolume.
C Example:
 status = RM_SetUnitsSolution(id, 1);
"""
function RM_SetUnitsSolution(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsSolution, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set input units for solid-solution assemblages. In PHREEQC, solid solutions are defined by moles of each component (Mp). RM_SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single SOLID_SOLUTION definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of a solid-solution component will be the same regardless of porosity. For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.

Parameters
id	The instance id returned from RM_Create.
option	Units option for solid solutions: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 status = RM_SetUnitsSSassemblage(id, 1);
"""
function RM_SetUnitsSSassemblage(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsSSassemblage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set input units for surfaces. In PHREEQC input, surfaces are defined by moles of surface sites (Mp). RM_SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume (RM_SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (RM_SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P).RV.

If a single SURFACE definition is used for cells with different initial porosity, the three options scale quite differently. For option 0, the number of moles of surface sites will be the same regardless of porosity. For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume. For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.

Parameters
id	The instance id returned from RM_Create.
option	Units option for surfaces: 0, 1, or 2.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_InitialPhreeqc2Module, RM_InitialPhreeqcCell2Module, RM_SetPorosity, RM_SetRepresentativeVolume.
C Example:
 status = RM_SetUnitsSurface(id, 1);
"""
function RM_SetUnitsSurface(id::Int, option::Int)
  IRM_RESULT=ccall((:RM_SetUnitsSurface, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, option)
end

"""
Set solution concentrations in the reaction cells based on the vector of aqueous species concentrations (species_conc). This method is intended for use with multicomponent-diffusion transport calculations, and RM_SetSpeciesSaveOn must be set to true. The list of aqueous species is determined by RM_FindComponents and includes all aqueous species that can be made from the set of components. The method determines the total concentration of a component by summing the molarities of the individual species times the stoichiometric coefficient of the element in each species. Solution compositions in the reaction cells are updated with these component concentrations.

Parameters
id	The instance id returned from RM_Create.
species_conc	Array of aqueous species concentrations. Dimension of the array is (nxyz, nspecies), where nxyz is the number of user grid cells (RM_GetGridCellCount), and nspecies is the number of aqueous species (RM_GetSpeciesCount). Concentrations are moles per liter.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_FindComponents, RM_GetSpeciesConcentrations, RM_GetSpeciesCount, RM_GetSpeciesD25, RM_GetSpeciesZ, RM_GetSpeciesName, RM_GetSpeciesSaveOn, RM_SetSpeciesSaveOn.
C Example:
 status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double*) malloc((size_t) (nxyz * nspecies * sizeof(double)));
...
status = RM_SpeciesConcentrations2Module(id, species_c);
status = RM_RunCells(id);
"""
function RM_SpeciesConcentrations2Module(id::Int, species_conc::Array{Float64,1})
  IRM_RESULT=ccall((:RM_SpeciesConcentrations2Module, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, species_conc)
end

"""
Determines the volume and density to use when converting from the reaction-module concentrations to transport concentrations (RM_GetConcentrations). Two options are available to convert concentration units: (1) the density and solution volume calculated by PHREEQC are used, or (2) the specified density (RM_SetDensity) and solution volume are defined by the product of saturation (RM_SetSaturation), porosity (RM_SetPorosity), and representative volume (RM_SetRepresentativeVolume). Transport models that consider density-dependent flow will probably use the PHREEQC-calculated density and solution volume (default), whereas transport models that assume constant-density flow will probably use specified values of density and solution volume. Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat. Density is only used when converting to transport units of mass fraction.

Parameters
id	The instance id returned from RM_Create.
tf	True indicates that the solution density and volume as calculated by PHREEQC will be used to calculate concentrations. False indicates that the solution density set by RM_SetDensity and the volume determined by the product of RM_SetSaturation, RM_SetPorosity, and RM_SetRepresentativeVolume, will be used to calculate concentrations retrieved by RM_GetConcentrations.
See also
RM_GetConcentrations, RM_SetDensity, RM_SetPorosity, RM_SetRepresentativeVolume, RM_SetSaturation.
C Example:
 status = RM_UseSolutionDensityVolume(id, 0);
"""
function RM_UseSolutionDensityVolume(id::Int, tf::Int)
  IRM_RESULT=ccall((:RM_UseSolutionDensityVolume, Lib_PhreeqcRM_path), Cint,
  (Cint,Cint), id, tf)
end

"""
Print a warning message to the screen and the log file.

Parameters
id	The instance id returned from RM_Create.
warn_str	String to be printed.
Return values
IRM_RESULT	0 is success, negative is failure (See RM_DecodeError).
See also
RM_OpenFiles, RM_LogMessage, RM_OutputMessage, RM_ScreenMessage, RM_ErrorMessage.
C Example:
 status = RM_WarningMessage(id, "Parameter is out of range, using default");
"""
function RM_WarningMessage(id::Int, warn_str::AbstractString)
  IRM_RESULT=ccall((:RM_WarningMessage, Lib_PhreeqcRM_path), Cint,
  (Cint,Cstring), id, warn_str)
end

export RM_Abort,
 RM_CloseFiles,
 RM_Concentrations2Utility,
 RM_Create,
 RM_CreateMapping,
 RM_DecodeError,
 RM_Destroy,
 RM_DumpModule,
 RM_ErrorMessage,
 RM_FindComponents,
 RM_GetBackwardMapping,
 RM_GetChemistryCellCount,
 RM_GetComponent,
 RM_GetComponentCount,
 RM_GetConcentrations,
 RM_GetDensity,
 RM_GetEndCell,
 RM_GetErrorString,
 RM_GetErrorStringLength,
 RM_GetFilePrefix,
 RM_GetGfw,
 RM_GetGridCellCount,
 RM_GetIPhreeqcId,
 RM_GetNthSelectedOutputUserNumber,
 RM_GetSaturation,
 RM_GetSelectedOutput,
 RM_GetSelectedOutputColumnCount,
 RM_GetSelectedOutputCount,
 RM_GetSelectedOutputHeading,
 RM_GetSelectedOutputRowCount,
 RM_GetSolutionVolume,
 RM_GetSpeciesConcentrations,
 RM_GetSpeciesCount,
 RM_GetSpeciesD25,
 RM_GetSpeciesName,
 RM_GetSpeciesSaveOn,
 RM_GetSpeciesZ,
 RM_GetStartCell,
 RM_GetThreadCount,
 RM_GetTime,
 RM_GetTimeConversion,
 RM_GetTimeStep,
 RM_InitialPhreeqc2Concentrations,
 RM_InitialPhreeqc2Module,
 RM_InitialPhreeqc2SpeciesConcentrations,
 RM_InitialPhreeqcCell2Module,
 RM_LoadDatabase,
 RM_LogMessage,
# RM_MpiWorker,
# RM_MpiWorkerBreak,
 RM_OpenFiles,
 RM_OutputMessage,
 RM_RunCells,
 RM_RunFile,
 RM_RunString,
 RM_ScreenMessage,
 RM_SetComponentH2O,
 RM_SetConcentrations,
 RM_SetCurrentSelectedOutputUserNumber,
 RM_SetDensity,
 RM_SetDumpFileName,
 RM_SetErrorHandlerMode,
 RM_SetFilePrefix,
# RM_SetMpiWorkerCallback,
# RM_SetMpiWorkerCallbackCookie,
 RM_SetPartitionUZSolids,
 RM_SetPorosity,
 RM_SetPressure,
 RM_SetPrintChemistryMask,
 RM_SetPrintChemistryOn,
 RM_SetRebalanceByCell,
 RM_SetRebalanceFraction,
 RM_SetRepresentativeVolume,
 RM_SetSaturation,
 RM_SetScreenOn,
 RM_SetSelectedOutputOn,
 RM_SetSpeciesSaveOn,
 RM_SetTemperature,
 RM_SetTime,
 RM_SetTimeConversion,
 RM_SetTimeStep,
 RM_SetUnitsExchange,
 RM_SetUnitsGasPhase,
 RM_SetUnitsKinetics,
 RM_SetUnitsPPassemblage,
 RM_SetUnitsSolution,
 RM_SetUnitsSSassemblage,
 RM_SetUnitsSurface,
 RM_SpeciesConcentrations2Module,
 RM_UseSolutionDensityVolume,
 RM_WarningMessage
end
