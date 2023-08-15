# DNA Unzipping Curve Calculator (ver.3)

## How to build the project

Create a folder for building the project:

>**mkdir build**  
>**cd build**

Build "BPEnergy.lib". This lib includes functions that calculate unzip energy of a DNA sequence.

>**g++ -c ../BPEnergy/BPEnergy.cpp -o BPEnergy.o**  
>**ar rcs BPEnergy.lib BPEnergy.o**

Build "UnzipLUT.lib". This lib contains LUTs which is the key of the program that boost the speed hundres of times. *J_SIZE* and *EXT_SIZE* directives controls the size of LUTs. The large the numbers, the higher the precision. The build of this lib can be extremely slow, dependends on the size of LUTs.

>**g++ -c ../UnzipLUT/UnzipLUT.cpp -std=c++20 -fconstexpr-ops-limit=100000000000 -DJ_SIZE=64 -DEXT_SIZE=64 -o UnzipLUT.o**  
>**ar rcs UnzipLUT.lib UnzipLUT.o**

Build "main.cpp" / link libraries.

>**g++ ../main.cpp -L. -lUnzipLUT -lBPEnergy -o UnzipDNA.exe**

Do some test using the example data:

>**UnzipDNA.exe ../exampleData.txt**

*CMAKE somehow keeps generating build tree for visual studio, but MSVC has issues generating constexpr...*  

The executable built can be run like "unzip.exe exampleData.txt", where "exampleData.txt" file contains the DNA sequence to be unzipped.  

I first wrote the program in Python, and realized that it was too slow for my application. To speed up the calculation, I use look-up table (LUT) to save information that supposes to be calculated in run time. These LUTs are saved in constexpr containers.  

It indeed took me some thinking to move majority of calculation overheads from run-time to compile time. After several attempts, I almost **"constexpred"** everyting that can be calculated ahead of time. As a result, this program is 100 times faster than the python program, 6 times faster than the Version 2 program. The only drawback is that the compile time is very long, thousands of times longer than the prototype c++ program (which does not calculate many constexpr variables).  

## DNA unzipping theory

For more information on DNA unzipping experiment and theory:

>[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS.
>
>[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal.
>
>[3] Huguet, Josep M., et al. (2010) PNAS.

