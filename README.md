<<<<<<< HEAD
# unzipDNA_v2
-build command:
**g++ -std=c++20 -o unzip.exe src/main.cpp src/Optimization.cpp**

-how to use:
**unzip.exe**
or
**unzip.exe exampleData.txt 2200 100 298 0.07**

the CMAKE somehow only generates build tree for visual studio, but MSVC has issues generating constexpr...
=======
# Calculate the DNA unzipping curve 

**build command (Windows)**

> g++ -std=c++20 -o a.exe main.cpp && a.exe

* CMAKE somehow only generates build tree for visual studio. I cannot set MAKE_*
*  but vs has some issue when generate constexpr container...*
>>>>>>> 051e763b0cc5a39a76db23bbd7b89d583fa6b4c8
