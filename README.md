# unzipDNA_v2
-build command:
**g++ -std=c++20 -o unzip.exe src/main.cpp src/Optimization.cpp**

-how to use:
**unzip.exe**
or
**unzip.exe exampleData.txt 2200 100 298 0.07**

the CMAKE somehow only generates build tree for visual studio, but MSVC has issues generating constexpr...
