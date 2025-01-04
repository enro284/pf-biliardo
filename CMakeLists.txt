cmake_minimum_required(VERSION 3.16)
project(pf-biliardo VERSION 0.2.0)

include_directories("src")

# abilita il supporto per i test, tra cui l'opzione BUILD_TESTING usata sotto
include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# richiedi l'uso di C++20, senza estensioni non-standard offerte dal compilatore usato 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# abilita warning
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Wimplicit-fallthrough -Wextra-semi -Wold-style-cast")

# abilita:
# - l'address sanitizer e l'undefined-behaviour sanitizer
# - _LIBCPP_DEBUG (utilizzato da clang)
# - _GLIBCXX_ASSERTIONS (utilizzato da gcc)
#in debug mode
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer -D_LIBCPP_DEBUG -D_GLIBCXX_ASSERTIONS")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer -D_LIBCPP_DEBUG -D_GLIBCXX_ASSERTIONS")

# LIBRARIES
find_package(SFML 2.5 COMPONENTS graphics REQUIRED)

add_library(mathematics src/mathematics.cpp)

add_library(kinematics src/kinematics.cpp)
target_link_libraries(kinematics mathematics)

add_library(statistics src/statistics.cpp)
target_link_libraries(statistics mathematics)

add_library(graphics src/graphics.cpp)
target_link_libraries(graphics kinematics)

# EXECUTABLES
add_executable(biliardo src/main.cpp)
target_link_libraries(biliardo kinematics statistics graphics sfml-graphics)

add_executable(multiple_particle_sim_csv src/main_csv.cpp)
target_link_libraries(multiple_particle_sim_csv kinematics)


# TESTS
# se il testing e' abilitato...
#   per disabilitare il testing, passare -DBUILD_TESTING=OFF a cmake durante la fase di configurazione
if (BUILD_TESTING)

  # aggiungi l'eseguibile mathematics.t
  add_executable(mathematics.t tests/mathematics.test.cpp)
  target_link_libraries(mathematics.t mathematics)
  # aggiungi l'eseguibile mathematics.t alla lista dei test
  add_test(NAME mathematics.t COMMAND mathematics.t)

  # aggiungi l'eseguibile kinematics.t
  add_executable(kinematics.t tests/kinematics.test.cpp)
  target_link_libraries(kinematics.t kinematics)
  # aggiungi l'eseguibile kinematics.t alla lista dei test
  add_test(NAME kinematics.t COMMAND kinematics.t)
  
  # aggiungi l'eseguibile statistics.t
  add_executable(statistics.t tests/statistics.test.cpp src/statistics.cpp)
  target_link_libraries(statistics.t statistics)
  # aggiungi l'eseguibile statistics.t alla lista dei test
  add_test(NAME statistics.t COMMAND statistics.t)

endif()

