add_executable(egh_rps
               mod_rps.f90
               mod_hmm.f90
               egh_rps.f90)
install(TARGETS egh_rps
        DESTINATION bin)
