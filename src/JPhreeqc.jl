module JPhreeqc

const Lib_PhreeqcRM_path=joinpath("/usr/local/lib", "libphreeqcrm.so")
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
  IRM_RESULT=ccall((:RM_CreateMapping, Lib_PhreeqcRM_path),
  Cint, (Cint, Ptr{Cdouble}), id, grid2chem)
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
function RM_Abort(id::Int, result::Int, err_str::String)
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
  IRM_RESULT=ccall((:RM_GetChemistryCellCount, Lib_PhreeqcRM_path), Cint,
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
function RM_ErrorMessage(id::Int, errstr::String)
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
  n_components=ccall((:RM_GetChemistryCellCount, Lib_PhreeqcRM_path), Cint,
  (Cint,), id)
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
function RM_GetBackwardMapping(id::Int, n::Int, list::Array{Int,1}, size::Array{Int,1})

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
function RM_GetComponent(id::Int, num::Int, chem_name::String, l::Int)
  IRM_RESULT=ccall((:RM_GetComponent, Lib_PhreeqcRM_path), Cint,
  (Cint, Cint, Cstring, Cint), id, num, chem_name, l)
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
function RM_GetEndCell(id::Int, ec::Array{Int,1})
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
function RM_GetErrorString(id::Int, errstr::String, l::Int)
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
  er_length=ccall((:RM_GetErrorString, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
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
function RM_GetFilePrefix(id::Int, prefix::String, l::Int)
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
function RM_GetGfw(id::Int, gfw::Array{Float6,1})
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
function RM_GetSelectedOutputHeading(id::Int, icol::Int, heading::String, length::Int)
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
function RM_GetSpeciesName(id::Int, i::Int, name::String, length::Int)
  =ccall((:RM_GetSpeciesName, Lib_PhreeqcRM_path), Cint,
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
  =ccall((:RM_GetSpeciesZ, Lib_PhreeqcRM_path), Cint,
  (Cint,Ptr{Cdouble}), id, z)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

"""

"""
function (id::Int)
  =ccall((:, Lib_PhreeqcRM_path), Cint,  (Cint,), id)
end

























function RM_Abort(id::Cint,result::Cint,err_str::Ptr{UInt8})
        ccall((:RM_Abort,RM_interface_C),IRM_RESULT,(Cint,Cint,Ptr{UInt8}),id,result,err_str)
end
function RM_CloseFiles(id::Cint)
        ccall((:RM_CloseFiles,RM_interface_C),IRM_RESULT,(Cint,),id)
    end
    function RM_Concentrations2Utility(id::Cint,c::Ptr{Cdouble},n::Cint,tc::Ptr{Cdouble},p_atm::Ptr{Cdouble})
        ccall((:RM_Concentrations2Utility,RM_interface_C),Cint,(Cint,Ptr{Cdouble},Cint,Ptr{Cdouble},Ptr{Cdouble}),id,c,n,tc,p_atm)
    end
    function RM_Create(nxyz::Cint,nthreads::Cint)
        ccall((:RM_Create,RM_interface_C),Cint,(Cint,Cint),nxyz,nthreads)
    end
    function RM_CreateMapping(id::Cint,grid2chem::Ptr{Cint})
        ccall((:RM_CreateMapping,RM_interface_C),IRM_RESULT,(Cint,Ptr{Cint}),id,grid2chem)
    end
    function RM_DecodeError(id::Cint,e::Cint)
        ccall((:RM_DecodeError,RM_interface_C),IRM_RESULT,(Cint,Cint),id,e)
    end
    function RM_Destroy(id::Cint)
        ccall((:RM_Destroy,RM_interface_C),IRM_RESULT,(Cint,),id)
    end
    function RM_DumpModule(id::Cint,dump_on::Cint,append::Cint)
        ccall((:RM_DumpModule,RM_interface_C),IRM_RESULT,(Cint,Cint,Cint),id,dump_on,append)
    end
    function RM_ErrorMessage(id::Cint,errstr::Ptr{UInt8})
        ccall((:RM_ErrorMessage,RM_interface_C),IRM_RESULT,(Cint,Ptr{UInt8}),id,errstr)
    end
    function RM_FindComponents(id::Cint)
        ccall((:RM_FindComponents,RM_interface_C),Cint,(Cint,),id)
    end)
    function RM_SetUnitsExchange(id::Cint,option::Cint)
        ccall((:RM_SetUnitsExchange,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsGasPhase(id::Cint,option::Cint)
        ccall((:RM_SetUnitsGasPhase,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsKinetics(id::Cint,option::Cint)
        ccall((:RM_SetUnitsKinetics,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsPPassemblage(id::Cint,option::Cint)
        ccall((:RM_SetUnitsPPassemblage,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsSolution(id::Cint,option::Cint)
        ccall((:RM_SetUnitsSolution,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsSSassemblage(id::Cint,option::Cint)
        ccall((:RM_SetUnitsSSassemblage,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SetUnitsSurface(id::Cint,option::Cint)
        ccall((:RM_SetUnitsSurface,RM_interface_C),IRM_RESULT,(Cint,Cint),id,option)
    end
    function RM_SpeciesConcentrations2Module(id::Cint,species_conc::Ptr{Cdouble})
        ccall((:RM_SpeciesConcentrations2Module,RM_interface_C),IRM_RESULT,(Cint,Ptr{Cdouble}),id,species_conc)
    end
    function RM_UseSolutionDensityVolume(id::Cint,tf::Cint)
        ccall((:RM_UseSolutionDensityVolume,RM_interface_C),IRM_RESULT,(Cint,Cint),id,tf)
    end
    function RM_WarningMessage(id::Cint,warn_str::Ptr{UInt8})
        ccall((:RM_WarningMessage,RM_interface_C),IRM_RESULT,(Cint,Ptr{UInt8}),id,warn_str)
    end

export RM_create, RM_CloseFiles, RM_Concentrations2Utility, RM_CreateMapping,
  RM_Abort,
end
