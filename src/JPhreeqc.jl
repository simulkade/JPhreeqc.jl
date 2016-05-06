module JPhreeqc

const Lib_PhreeqcRM_path=joinpath("/usr/local/lib", "libphreeqcrm.so")
#const phreeqc_path1=Libdl.find_library(["libphreeqcrm"], ["/usr/local/lib"])

#n_cell=1
#n_thread=1
function RM_Create(n_cell::Int, n_thread::Int)
  phreeqc_id= ccall((:RM_Create, Lib_PhreeqcRM_path), Cint, (Cint, Cint), n_cell, n_thread)
end

export RM_create
end
