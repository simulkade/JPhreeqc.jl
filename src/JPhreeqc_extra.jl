# some extra utility functions for an easier use of PhreeqcRM in Julia

"""
components_list=getComponentList(id)
returns the list of components (elements and charge) as a cell array of strings
"""
function getComponentList(id::Int)
  ncomps = RM_FindComponents(id)
  components = cell(ncomps)
	for i = 1:ncomps
		components[i] = string(zeros(10))
		status = RM_GetComponent(id, i-1, components[i], length(components[i]))
    components[i]=strip(components[i], '\0')
	end
  return components
end

"""
species_list=getSpeciesList(id)
returns the list of ionic species as a cell array of strings
"""
function getSpeciesList(id::Int)
  ncomps = RM_FindComponents(id)
  n_species = RM_GetSpeciesCount(id)
  species = cell(n_species)
	for i = 1:n_species
		species[i] = string(zeros(10))
		status = RM_GetComponent(id, i-1, species[i], length(species[i]))
    species[i]=strip(species[i], '\0')
	end
  return species
end

function setDefaultPhreeqcProperties(id::Int)
  # Set properties
  status = RM_SetErrorHandlerMode(id, 2)
	status = RM_SetComponentH2O(id, 0)
	status = RM_SetRebalanceFraction(id, 0.5)
	status = RM_SetRebalanceByCell(id, 1)
	status = RM_UseSolutionDensityVolume(id, 0)
	status = RM_SetPartitionUZSolids(id, 0)
end

function setDefaultPhreeqcUnits(id::Int)
  status = RM_SetUnitsSolution(id, 2)      # 1, mg/L 2, mol/L 3, kg/kgs
	status = RM_SetUnitsPPassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsExchange(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSurface(id, 1)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsGasPhase(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSSassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsKinetics(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
end

"""
The default call to this function assumes that an aqueous solution exists in all cells.
aq_solution: 1 means exist, -1 means does not exist
"""
function setInitialVectors(id::Int, nxyz::Int; aq_solution::Int=1, eq_phase::Int=0, ion_exchange::Int=0,
                             surface_site::Int==0, gas_phase::Int=0, solid_solution::Int=0, kin_reaction::Int=0) 
	ic1 = zeros(Int, nxyz, 7)
	ic2 = zeros(Int, nxyz, 7)
	f1 = zeros(Float64, nxyz, 7)
	ic1[:,1] = 1       # Solution 1
	ic1[:,2] = -1      # Equilibrium phases none
	ic1[:,3] = -1       # Exchange none
	ic1[:,4] = -1      # Surface 1
	ic1[:,5] = -1      # Gas phase none
	ic1[:,6] = -1      # Solid solutions none
	ic1[:,7] = -1      # Kinetics none
	ic2[:,1] = -1      # Solution none
	ic2[:,2] = -1      # Equilibrium phases none
	ic2[:,3] = -1      # Exchange none
	ic2[:,4] = -1      # Surface none
	ic2[:,5] = -1      # Gas phase none
	ic2[:,6] = -1      # Solid solutions none
	ic2[:,7] = -1      # Kinetics none
	f1[:,1]  = 1.0      # Mixing fraction ic1 Solution
	f1[:,2] = 1.0      # Mixing fraction ic1 Equilibrium phases
	f1[:,3] = 1.0      # Mixing fraction ic1 Exchange 1
	f1[:,4] = 1.0      # Mixing fraction ic1 Surface
	f1[:,5] = 1.0      # Mixing fraction ic1 Gas phase
	f1[:,6] = 1.0      # Mixing fraction ic1 Solid solutions
	f1[:,7] = 1.0      # Mixing fraction ic1 Kinetics
	return ic1, ic2, f1
end

"""
id: phreeqcRM instance id
isel: selected output block number
nxyx: number of grids
"""
function getSelectedOutputArray(id::Int, nxyz::Int, isel::Int)
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
  status = RM_SetCurrentSelectedOutputUserNumber(id, Int(n_user))
  # Get double array of selected output values
  col = RM_GetSelectedOutputColumnCount(id)
  # allocate(selected_out(nxyz,col))
  selected_out = zeros(Float64, nxyz*col)
  status = RM_GetSelectedOutput(id, selected_out)
  return reshape(selected_out, nxyz, col)
end

function getSelectedOutputHeading(id::Int, isel::Int)
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
  status = RM_SetCurrentSelectedOutputUserNumber(id, Int(n_user))
  col = RM_GetSelectedOutputColumnCount(id)
  heading=cell(col)
  for j = 0:col-1
    heading[j+1]=string(zeros(Int, 25))
    status = RM_GetSelectedOutputHeading(id, j, heading[j+1], length(heading[j+1]))
    heading[j+1]=strip(heading[j+1], '\0')
  end
  return heading
end

