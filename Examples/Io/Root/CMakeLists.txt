# add_definitions(-g -O2)
add_library(
  ActsExamplesIoRoot SHARED
  src/RootMeasurementWriter.cpp
  src/RootMaterialDecorator.cpp
  src/RootMaterialWriter.cpp
  src/RootMaterialTrackReader.cpp
  src/RootMaterialTrackWriter.cpp
  src/RootPlanarClusterWriter.cpp
  src/RootParticleWriter.cpp
  src/RootParticleReader.cpp
  src/RootPropagationStepsWriter.cpp
  src/RootSimHitWriter.cpp
  src/RootSpacepointWriter.cpp
  src/RootTrackParameterWriter.cpp
  src/RootTrajectoryStatesWriter.cpp
  src/RootTrajectorySummaryReader.cpp
  src/RootTrajectorySummaryWriter.cpp
  src/RootVertexPerformanceWriter.cpp
  src/RootBFieldWriter.cpp
  src/RootAthenaNTupleReader.cpp
  src/RootDigiBase.cpp
)

target_include_directories(
  ActsExamplesIoRoot
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoRoot
  PUBLIC
    ActsCore ActsPluginIdentification ActsExamplesDigitization
    ActsExamplesFramework ActsExamplesPropagation ActsExamplesTruthTracking ActsExamplesMaterialMapping
    Threads::Threads
  PRIVATE ROOT::Core ROOT::Hist ROOT::Tree)

install(
  TARGETS ActsExamplesIoRoot
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

root_generate_dictionary(ActsExamplesIoRootDict  MODULE ActsExamplesIoRoot LINKDEF src/ActsExamplesIoRootLinkDef.h)
