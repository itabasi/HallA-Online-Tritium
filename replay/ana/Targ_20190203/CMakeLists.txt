#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(HRSL_nnL)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
#set(ROOT_INCLUDE_DIR /home/gogami/root/root-5.34.36-install/include/root)
#set(ROOT_INCLUDE_DIR /usr/local/include)
#set(ROOT_LIBRARIES -L/usr/local/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic)

############## Itabashi ##############
set(ROOT_INCLUDE_DIR /usr/local/include)
set(ROOT_LIBRARIES -L/usr/local/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
add_executable(TARG TARG.cc ${sources} ${headers})
target_link_libraries(TARG ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
#set(CMAKE_INSTALL_PREFIX /home/gogami/11bgogami/svnroot_JLab/E12-17-003/simulation/TARG)
set(CMAKE_INSTALL_PREFIX /home/dragon/HD-PZU3/dragon/svnroot_JLab/E12-17-003/simulations/TARG/)
install(TARGETS TARG RUNTIME DESTINATION bin)



