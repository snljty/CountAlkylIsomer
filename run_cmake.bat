@echo off

set CMAKE_PREFIX_PATH=
cmake . -B build -G "MinGW Makefiles" -D CMAKE_INSTALL_PREFIX=%CD% -LH
rem cmake --build build -j
rem cmake --install build
