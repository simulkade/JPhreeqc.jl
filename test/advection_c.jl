"""
advection_c example from PhreeqcRM test cases, converted to Julia
advection is solved using JFVM package, hence some changes in mapping cells
phreeqc is called via JPhreeqc package
written by A.A. Eftekhari, 8 May 2016
An exceptionally beautiful and sunny afternoon
Lyngby, Denmark
before using this function, write:
using JFVM, JPhreeqc
"""
function advection_c()
  Nx=40
  Ny=1
  Nz=1
  nxyz=Nx*Ny*Nz # number of cells
  nthreads=1 # number of threads
  id=RM_Create(nxyz, nthreads) # create a phreeqc instance
  # Set properties
  status = RM_SetErrorHandlerMode(id, 2)
	status = RM_SetComponentH2O(id, 0)
	status = RM_SetRebalanceFraction(id, 0.5)
	status = RM_SetRebalanceByCell(id, 1)
	status = RM_UseSolutionDensityVolume(id, 0)
	status = RM_SetPartitionUZSolids(id, 0)
	status = RM_SetFilePrefix(id, "Advect_c")
  status = RM_OpenFiles(id)
  # Set concentration units
  status = RM_SetUnitsSolution(id, 2)      # 1, mg/L 2, mol/L 3, kg/kgs
	status = RM_SetUnitsPPassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsExchange(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSurface(id, 1)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsGasPhase(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsSSassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
	status = RM_SetUnitsKinetics(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
  # Set conversion from seconds to days
  # I prefer to work with seconds, so I keep it to 1.0
  status = RM_SetTimeConversion(id, 1.0)
  # Set representative volume
  rv = ones(Float64, nxyz) # volume is in liter
  status = RM_SetRepresentativeVolume(id, rv)
  # Set initial porosity
  por = 1.0*ones(Float64, nxyz)
	status = RM_SetPorosity(id, por)
  # Set initial saturation
	sat = ones(Float64, nxyz)
	status = RM_SetSaturation(id, sat)
	# Set cells to print chemistry when print chemistry is turned on
	print_chemistry_mask = ones(Int32, nxyz)
	status = RM_SetPrintChemistryMask(id, print_chemistry_mask)
	# Partitioning of uz solids
	status = RM_SetPartitionUZSolids(id, 0)
	# Demonstation of mapping, no symmetry, only one row of cells
	grid2chem = collect(Int32, 0:nxyz)
	status = RM_CreateMapping(id, grid2chem)
	if (status < 0)
    status = RM_DecodeError(id, status)
  end
	nchem = RM_GetChemistryCellCount(id)

	# --------------------------------------------------------------------------
	# Set initial conditions
	# --------------------------------------------------------------------------
	# Set printing of chemistry file
	status = RM_SetPrintChemistryOn(id, 0, 1, 0) # workers, initial_phreeqc, utility
	# Set printing of chemistry file
	status = RM_LoadDatabase(id, "/home/ali/MyPackages/JPhreeqc.jl/test/phreeqc.dat")

	# Demonstration of error handling if ErrorHandlerMode is 0
  # commented out for now (AAE)
	# if (status != IRM_OK)
	# {
	# 	l = RM_GetErrorStringLength(id)
	# 	errstr = (char *) malloc((size_t) (l * sizeof(char) + 1))
	# 	RM_GetErrorString(id, errstr, l+1)
	# 	fprintf(stderr,"Beginning of error string:\n")
	# 	fprintf(stderr,"%s", errstr)
	# 	fprintf(stderr,"End of error string.\n")
	# 	free(errstr)
	# 	errstr = NULL
	# 	RM_Destroy(id)
	# 	exit(1)
	# }
	# Run file to define solutions and reactants for initial conditions, selected output
	# There are three types of IPhreeqc instances in PhreeqcRM
	# Argument 1 refers to the workers for doing reaction calculations for transport
	# Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
	# Argument 3 refers to the Utility instance
	status = RM_RunFile(id, 1, 1, 1, "/home/ali/MyPackages/JPhreeqc.jl/test/advect.pqi")
	# Clear contents of workers and utility
	str="DELETE -all"
	status = RM_RunString(id, 1, 0, 1, str)	# workers, initial_phreeqc, utility
	# Determine number of components to transport
	ncomps = RM_FindComponents(id)
	# Print some of the reaction module information
  # commented out by AAE
	# sprintf(str1, "Number of threads:                                %d\n", RM_GetThreadCount(id))
	# status = RM_OutputMessage(id, str1)
	# sprintf(str1, "Number of MPI processes:                          %d\n", RM_GetMpiTasks(id))
	# status = RM_OutputMessage(id, str1)
	# sprintf(str1, "MPI task number:                                  %d\n", RM_GetMpiMyself(id))
	# status = RM_OutputMessage(id, str1)
	# status = RM_GetFilePrefix(id, str, 100)
	# sprintf(str1, "File prefix:                                      %s\n", str)
	# status = RM_OutputMessage(id, str1)
	# sprintf(str1, "Number of grid cells in the user's model:         %d\n", RM_GetGridCellCount(id))
	# status = RM_OutputMessage(id, str1)
	# sprintf(str1, "Number of chemistry cells in the reaction module: %d\n", RM_GetChemistryCellCount(id))
	# status = RM_OutputMessage(id, str1)
	# sprintf(str1, "Number of components for transport:               %d\n", RM_GetComponentCount(id))
	# status = RM_OutputMessage(id, str1)
	# Get component information
	components = cell(ncomps)
	gfw = zeros(Float64, ncomps)
	status = RM_GetGfw(id, gfw)
	for i = 1:ncomps
		components[i] = string(zeros(10))
		status = RM_GetComponent(id, i-1, components[i], length(components[i]))
		print(components[i], "\t", gfw[i], "\n")
	end
	# Set array of initial conditions
  # these arrays must be 1D, however it is easier to define them as 2D
  # and use the `:` operator to convert it to 1D
  # each column is for different phreeqc phases
	ic1 = zeros(Int32, nxyz, 7)
	ic2 = zeros(Int32, nxyz, 7)
	f1 = zeros(Float64, nxyz, 7)
	ic1[:,1] = 1       # Solution 1
	ic1[:,2] = -1      # Equilibrium phases none
	ic1[:,3] = 1       # Exchange 1
	ic1[:,4] = -1      # Surface none
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
	status = RM_InitialPhreeqc2Module(id, ic1[:], ic2[:], f1[:])
	# No mixing is defined, so the following is equivalent
  # AAE: the following does not work with the current JPhreeqc
	# status = RM_InitialPhreeqc2Module(id, ic1, NULL, NULL)

	# alternative for setting initial conditions
	# cell number in second argument (-1 indicates last solution, 40 in this case)
	# in advect.pqi and any reactants with the same number--
	# Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
	# will be written to cells 18 and 19 (0 based)
	#module_cells = (int *) malloc((size_t) (2 * sizeof(int)))
	#module_cells[0] = 18
	#module_cells[1] = 19
	#status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2)
  # AAE: I could not really understant what happens above!
	# Initial equilibration of cells
	time1 = 0.0 # changed time to time1 to avoid conflicts with julia time function
	time_step = 0.0
	c = zeros(Float64, ncomps*nxyz)
	status = RM_SetTime(id, time1)
	status = RM_SetTimeStep(id, time_step)
	status = RM_RunCells(id)
	status = RM_GetConcentrations(id, c) # this is the most important line!
  # AAE: what I did in matlab was to get the volume from phreeqc and scale the
  # amount of solution in each cell
  # print(c)
	# --------------------------------------------------------------------------
	# Set boundary condition
	# --------------------------------------------------------------------------

	nbound = 1
	bc1 = zeros(Int32,nbound)
  bc2 = zeros(Int32,nbound)
	bc_f1 = zeros(Float64, nbound)
	bc_conc = zeros(Float64, ncomps*nbound)
	bc1[:]          = 0       # Solution 0 from Initial IPhreeqc instance
	bc2[:]          = -1      # no bc2 solution for mixing
	bc_f1[:]        = 1.0     # mixing fraction for bc1

	status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1)
  print(bc_conc)
	# --------------------------------------------------------------------------
	# Transient loop
	# --------------------------------------------------------------------------
  # this is the beautiful part: I can easily set values after the transport step

	nsteps = 10
	density = zeros(Float64, nxyz)
	volume = zeros(Float64, nxyz)
	pressure = zeros(Float64, nxyz)
	temperature = zeros(Float64, nxyz)
	sat_calc = zeros(Float64, nxyz)
	density[:] = 1.0
	pressure[:] = 2.0
	temperature[:] = 20.0

	status = RM_SetDensity(id, density)
	status = RM_SetPressure(id, pressure)
	status = RM_SetTemperature(id, temperature)
	time_step = 720.0
	status = RM_SetTimeStep(id, time_step)
	# for (isteps = 0 isteps < nsteps isteps++)
	# {
	# 	# Advection calculation
	# 	sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ",
	# 		RM_GetTime(id) * RM_GetTimeConversion(id), " days\n")
	# 	status = RM_LogMessage(id, str)
	# 	status = RM_SetScreenOn(id, 1)
	# 	status = RM_ScreenMessage(id, str)
	# 	sprintf(str, "%s%10.1f%s", "          Time step                  ",
	# 		RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days\n")
	# 	status = RM_LogMessage(id, str)
	# 	status = RM_ScreenMessage(id, str)

# ============================================================================
    # My JFVM code goes here!
    n_cells=nxyz
    cell_length=0.002 # [m]
    n_shifts=100*n_cells/40
    ncomp=ncomps
    c_effluent=zeros(125,ncomp)
    time_step=720.0 # [s]
    t_final=n_shifts*time_step # [s]
    porosity=1.0 # only for this example
    u_x = cell_length/time_step*porosity
    Lx=n_cells*cell_length # [m] domain length
    m=createMesh1D(Nx, Lx)
    Ux=createFaceVariable(m, u_x)
    BC=cell(ncomp)
    BC_right=zeros(ncomp)
    for i=1:ncomp
      BC[i]=createBC(m)
      # Left: Dirichlet
      BC[i].left.a[:]=0.0
      BC[i].left.b[:]=1.0
      BC[i].left.c[:]=bc_conc[i] # comes from phreeqc
      # Right: Neumann
      BC[i].right.a[:]=1.0
      BC[i].right.b[:]=0.0
      BC[i].right.c[:]=BC_right[i]
    end
    # n_steps=20
    # t_end=Lx/u_x
    # dt=t_end/20
    # initial conditions
    c_init=reshape(c,n_cells, Int64(ncomp))
    c_transport=zeros(Float64, n_cells, ncomp)
    C_old=cell(ncomp)
    C_new=cell(ncomp)
    for i=1:ncomp
      C_old[i]=createCellVariable(m, c_init[:,i], BC[i]) # only for 1D; reshapeCell for 2D,3D
      C_new[i]=copyCell(C_old[i])
    end
    M_bc     = cell(ncomp)
    RHS_bc   = cell(ncomp)
    M_conv   = cell(ncomp)
    RHS_conv = cell(ncomp)
    M_t      = cell(ncomp)
    RHS_t    = cell(ncomp)
    for i=1:ncomp
        M_bc[i], RHS_bc[i]=boundaryConditionTerm(BC[i])
        M_conv[i]=convectionUpwindTerm(Ux)
    end

    dt=time_step
    t_end=t_final
    t_step=0
    for t=dt:dt:t_end
      for i=1:ncomp
        M_t[i], RHS_t[i]=transientTerm(C_old[i], dt)
      end
      for i=1:ncomp
          C_new[i]=solveLinearPDE(m, M_t[i]+M_conv[i]+M_bc[i],
                  RHS_t[i]+RHS_bc[i]) # add +RHS_conv{i} to make it TVD
          c_transport[:,i]=C_new[i].value[2:end-1]
      end
    #advect_c(c, bc_conc, ncomps, nxyz, nbound)
    # ============================================================================

      #AAE: nothing changes in this code, so commented out
  		# Transfer data to PhreeqcRM for reactions
  		# status = RM_SetPorosity(id, por)              # If porosity changes
  		# status = RM_SetSaturation(id, sat)            # If saturation changes
  		# status = RM_SetTemperature(id, temperature)   # If temperature changes
  		# status = RM_SetPressure(id, pressure)         # If pressure changes
  		status = RM_SetConcentrations(id, c_transport[:])          # Transported concentrations
  	  status = RM_SetTimeStep(id, time_step)        # Time step for kinetic reactions
  		time1 = time1+time_step
  		status = RM_SetTime(id, time1)                 # Current time

  		# Run cells with transported conditions
  		# sprintf(str, "%s%10.1f%s", "Beginning reaction calculation       ", RM_GetTime(id) * RM_GetTimeConversion(id), " days\n")
  		# status = RM_LogMessage(id, str)
  		# status = RM_ScreenMessage(id, str)
      #print("reached this point")
      #sleep(20)
  		status = RM_RunCells(id)
  		# Transfer data from PhreeqcRM for transport
  		status = RM_GetConcentrations(id, c)          # Concentrations after reaction
      c_reac=reshape(c,n_cells, Int64(ncomp))
      print(c_transport[end-1,4], "\t", c_reac[end-1,4], "\n")
      t_step=t_step+1
      c_effluent[t_step,:]=c_reac[end,:]
      for i=1:ncomp
        C_old[i]=createCellVariable(m, c_reac[:,i], BC[i]) # only for 1D; reshapeCell for 2D,3D
      end
    end

  # Finalize
  status = RM_CloseFiles(id)
  #status = RM_MpiWorkerBreak(id)
  status = RM_Destroy(id)
  c_effluent
  # t_step
#end # function sdvection_c
end
