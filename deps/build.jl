if Sys.iswindows()
    url = "https://github.com/simulkade/PhreeqcRM/raw/master/lib/PhreeqcRM.dll"
    try
        println("Downloading the library files ...")
        download(url, joinpath(@__DIR__, "..", "lib", "PhreeqcRM.dll"))
    catch
        @warn "Could not downlaod the library file PhreeqcRM.dll"
        @warn "Down load it from: https://github.com/simulkade/PhreeqcRM/raw/master/lib/PhreeqcRM.dll"
        @warn "and copy it into the JPhreeqc/lib folder"
    end
elseif Sys.islinux()
    try
        url = "https://github.com/simulkade/PhreeqcRM/raw/master/linux_lib/libphreeqcrm-3.5.0.so"
        download(url, joinpath(@__DIR__, "..", "lib", "libphreeqcrm.so"))
    catch
        @warn "Could not downlaod the library file libphreeqcrm.so"
        @warn "Down load it from: https://github.com/simulkade/PhreeqcRM/raw/master/linux_lib/libphreeqcrm-3.5.0.so"
        @warn "and copy it into the JPhreeqc/lib folder"
    end
else
    @error "No binary files for your OS. Only Windows and Linux at the moment!"
end
