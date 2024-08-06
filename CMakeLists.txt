cmake_minimum_required(VERSION 3.16)
project(pf-biliardo VERSION 0.1.0)

# abilita il supporto per i test, tra cui l'opzione BUILD_TESTING usata sotto
include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# richiedi l'uso di C++17, senza estensioni non-standard offerte dal compilatore usato 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# abilita warning
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

# abilita:
# - l'address sanitizer e l'undefined-behaviour sanitizer
# - _LIBCPP_DEBUG (utilizzato da clang)
# - _GLIBCXX_ASSERTIONS (utilizzato da gcc)
#in debug mode
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer -D_LIBCPP_DEBUG -D_GLIBCXX_ASSERTIONS")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer -D_LIBCPP_DEBUG -D_GLIBCXX_ASSERTIONS")

add_executable(single_particle_sim main_one.cpp mathematics.cpp kinematics.cpp)
#add_executable(multiple_particle_sim main_mult.cpp mathematics.cpp kinematics.cpp)

# se il testing e' abilitato...
#   per disabilitare il testing, passare -DBUILD_TESTING=OFF a cmake durante la fase di configurazione
if (BUILD_TESTING)

  # aggiungi l'eseguibile mathematics.t
  add_executable(mathematics.t mathematics.test.cpp mathematics.cpp)
  # aggiungi l'eseguibile mathematics.t alla lista dei test
  add_test(NAME mathematics.t COMMAND mathematics.t)
  # aggiungi l'eseguibile kinematics.t
  add_executable(kinematics.t kinematics.test.cpp kinematics.cpp mathematics.cpp)
  # aggiungi l'eseguibile kinematics.t alla lista dei test
  add_test(NAME kinematics.t COMMAND kinematics.t)
  # aggiungi l'eseguibile statistics.t
  add_executable(statistics.t statistics.test.cpp statistics.cpp)
  # aggiungi l'eseguibile statistics.t alla lista dei test
  add_test(NAME statistics.t COMMAND statistics.t)

endif()

