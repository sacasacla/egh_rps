cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE 1)
enable_testing()
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  message(STATUS "Using GNU Fortran compiler")
  set(dialect "-O3 -Wuninitialized -ffree-form -std=f2008 -fimplicit-none")
  set(bounds "-fbounds-check -fbacktrace")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Using intel Fortran compiler")
  set(dialect "-O0 -free -implicitnone -qopenmp -mkl -heap-arrays -qopenmp")
  set(bounds "-traceback -check bounds -warn all")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

add_subdirectory(src)

add_test( NAME test_isonoke_02_05
  COMMAND egh_rps ${CMAKE_CURRENT_SOURCE_DIR}/dat/isonoke.txt 2 5
  )

add_test( NAME test_isonoke_10_05
  COMMAND egh_rps ${CMAKE_CURRENT_SOURCE_DIR}/dat/isonoke.txt 10 5
  )
