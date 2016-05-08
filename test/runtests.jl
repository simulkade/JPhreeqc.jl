using JPhreeqc
using Base.Test

# write your own tests here
phqc_id=RM_Create(1,1)
RM_LoadDatabase(phqc_id, "phreeqc_oil.dat")
RM_RunFile(phqc_id, 1,1,1,"chalkoil.pqi")
RM_Destroy(phqc_id)
@test 1==1
