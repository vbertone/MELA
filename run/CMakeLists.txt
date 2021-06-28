# build tests
file(GLOB testcodes *.f)
foreach(testsource ${testcodes})
  GET_FILENAME_COMPONENT(filename ${testsource} NAME_WE)
  add_executable(${filename} ${filename}.f)
  target_link_libraries(${filename} mela)
  add_test(${filename} ${filename})
endforeach()
