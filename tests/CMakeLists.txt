cmake_minimum_required (VERSION 3.9)

project (dephier)

add_executable(richdem_tests.exe
  tests.cpp
)

target_link_libraries(richdem_tests.exe
  PRIVATE 
    richdem
)
