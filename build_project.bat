@ echo off
set /P LUT_SIZE="Enter look up table size (Default=128): "
if "%LUT_SIZE%" == "" set "LUT_SIZE=128"

cd /d "%~dp0"
if not exist build mkdir build
cd build

for %%f in (*.o *.lib) do (
    del "%%f"
)

set start_time=%TIME%

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
g++ ../main.cpp -L. -lUnzipLUT -lBPEnergy -o UnzipDNA_%LUT_SIZE%x%LUT_SIZE%.exe

set end_time=%TIME%
for /f "tokens=1-3 delims=:." %%a in ("%start_time%") do (
    set /a "start_seconds=(((%%a*60)+%%b)*60)+%%c"
)
for /f "tokens=1-3 delims=:." %%a in ("%end_time%") do (
    set /a "end_seconds=(((%%a*60)+%%b)*60)+%%c"
)

REM Calculate time difference in seconds
set /a "time_difference=end_seconds-start_seconds"
set /a "hours=time_difference/3600"
set /a "minutes=(time_difference%%3600)/60"
set /a "seconds=time_difference%%60"
set "time_difference=%hours% hour %minutes% minute %seconds% second."
echo Compiled in: %time_difference% 

echo Testing the program..
UnzipDNA_%LUT_SIZE%x%LUT_SIZE%.exe ../exampleData.txt ..//out_exampleData_%LUT_SIZE%x%LUT_SIZE%.csv

cd ..