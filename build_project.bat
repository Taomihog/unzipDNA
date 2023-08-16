@ echo off
set /P LUT_SIZE="Enter look up table size (Default=128): "
if "%LUT_SIZE%" == "" set "LUT_SIZE=128"

cd /d "%~dp0"
if not exist build mkdir build
cd build

echo Building BPEnergy.lib...
REM basepair energy related stuff
g++ -c ../BPEnergy/BPEnergy.cpp -o BPEnergy.o
ar rcs BPEnergy.lib BPEnergy.o

echo Building UnzipLUT.lib...
REM create LUT
g++ -c ../UnzipLUT/UnzipLUT.cpp -std=c++20 -fconstexpr-ops-limit=100000000000 -DJ_SIZE=%LUT_SIZE% -DEXT_SIZE=%LUT_SIZE% -o UnzipLUT.o
ar rcs UnzipLUT.lib UnzipLUT.o

echo Building main...
REM executable
g++ ../main.cpp -L. -lUnzipLUT -lBPEnergy -o UnzipDNA.exe

echo Testing the program..
UnzipDNA.exe ../exampleData.txt

cd ..