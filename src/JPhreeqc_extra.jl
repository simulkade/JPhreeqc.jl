# some extra utility functions for an easier use of PhreeqcRM in Julia

"""
components_list=getComponentList(id)
returns the list of components as a cell array of strings
"""
function getComponentList(id)
  ncomps = RM_FindComponents(id)
  components = cell(ncomps)
	for i = 1:ncomps
		components[i] = string(zeros(10))
		status = RM_GetComponent(id, i-1, components[i], length(components[i]))
	end
  return components
end

function setDefaultPhreeqcProperties(id)
  # Set properties
  status = RM_SetErrorHandlerMode(id, 2)
	status = RM_SetComponentH2O(id, 0)
	status = RM_SetRebalanceFraction(id, 0.5)
	status = RM_SetRebalanceByCell(id, 1)
	status = RM_UseSolutionDensityVolume(id, 0)
	status = RM_SetPartitionUZSolids(id, 0)
end

function setDefaultPhreeqcUnits(id)
  status = RM_SetUnitsSolution(id, 2)      # 1, mg/L 2, mol/L 3, kg/kgs
	status = RM_SetUnitsPPassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsExchange(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSurface(id, 1)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsGasPhase(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSSassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsKinetics(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
end
