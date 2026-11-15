add_executable(spmm domain_melt5.f90)
add_executable(spmm_bond spmm_bond.f90)   # surface-energy prototype, see README
if(SPMM_NATIVE AND native)
  target_compile_options(spmm PRIVATE ${native})
  target_compile_options(spmm_bond PRIVATE ${native})
endif()
install(TARGETS spmm spmm_bond RUNTIME DESTINATION bin)
