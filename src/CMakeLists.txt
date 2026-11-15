add_executable(spmm spmm.f90)
if(SPMM_NATIVE AND native)
  target_compile_options(spmm PRIVATE ${native})
endif()
install(TARGETS spmm RUNTIME DESTINATION bin)

# Deprecated earlier program (centre-of-mass movement rule), not built by default.
if(SPMM_BUILD_LEGACY)
  add_executable(spmm_legacy domain_melt5.f90)
  if(SPMM_NATIVE AND native)
    target_compile_options(spmm_legacy PRIVATE ${native})
  endif()
  install(TARGETS spmm_legacy RUNTIME DESTINATION bin)
endif()
