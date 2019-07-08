using JPhreeqc
using Base.Test
include("test_phreeqc_adv.jl")
# write your own tests here
phqc_id=RM_Create(1,1)
RM_LoadDatabase(phqc_id, "phreeqc_oil.dat")
RM_RunFile(phqc_id, 1,1,1,"chalkoil.pqi")
RM_Destroy(phqc_id)
test_advect()
@test 1==1
