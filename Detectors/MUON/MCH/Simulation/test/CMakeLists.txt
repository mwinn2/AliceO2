O2_SETUP(NAME MCHSimulationTest)
set(BUCKET_NAME mch_simulation_test_bucket)


O2_GENERATE_EXECUTABLE(
  EXE_NAME test_MCHGeometry
  SOURCES testGeometry.cxx
  BUCKET_NAME mch_simulation_test_bucket
  )

O2_GENERATE_EXECUTABLE(
  EXE_NAME test_MCHSimulationResponse
  SOURCES testResponse.cxx
  BUCKET_NAME mch_simulation_test_bucket
  )

target_link_libraries(test_MCHGeometry Boost::unit_test_framework MCHSimulation)
add_test(NAME test_MCHGeometry COMMAND test_MCHGeometry)
target_link_libraries(test_MCHSimulationResponse Boost::unit_test_framework MCHSimulation)
add_test(NAME test_MCHSimulationResponse COMMAND test_MCHSimulationResponse)
