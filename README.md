# DNA Unzipping Curve Calculator (ver.3)


*CMAKE somehow keeps generating build tree for visual studio, but MSVC has issues generating constexpr...*   
I had to build by command line. 

## How to build the project

Create a folder for build:

>mkdir build
>
>cd build

Build "BPEnergy.lib". This lib includes functions that calculate unzip energy of a DNA sequence.


>g++ -c ../BPEnergy/BPEnergy.cpp -o BPEnergy.o
>
>ar rcs BPEnergy.lib BPEnergy.o

  
Build "UnzipLUT.lib". This lib contains LUTs which reduces the run-time calculation overhead and boosts the execution speed by a hundred of times.  

**J_SIZE** and **EXT_SIZE** control the sizes of LUTs. The large the numbers, the higher the precision.  

The build of this lib can be extremely slow, dependends on the size of LUTs.


>g++ -c ../UnzipLUT/UnzipLUT.cpp -std=c++20 -fconstexpr-ops-limit=100000000000 -DJ_SIZE=64 -DEXT_SIZE=64 -o UnzipLUT.o
>
>ar rcs UnzipLUT.lib UnzipLUT.o

  
Build "main.cpp" / link libraries.


>g++ ../main.cpp -L. -lUnzipLUT -lBPEnergy -o UnzipDNA.exe


Do some test using the example data:


>UnzipDNA.exe ../exampleData.txt  


The executable built can be run like "unzip.exe exampleData.txt", where "exampleData.txt" file contains the DNA sequence to be unzipped.  
  
  
It took some thinking to move majority of calculation from run-time to compile time. After several attempts, I **"constexpred"** most of the calculation overhead.   
  
Two look-up tables (LUTs) are used to hold these data. These LUTs are saved in **constexpr std::arrays**. The drawback is that the compile time is very long, thousands of times longer than the prototype c++ program (which does not calculate many constexpr variables).  
  
On Aug/15/2023, I implemented multithreading. The execution speed increased by another factor of 10-20.  
  
*At the current stage, this program is >1000 times faster than my original python program.*
  
## DNA unzipping theory
  
**DNA unzipping experiment on a 4.4 kb DNA**. single-molecule measurement (blue) and computer-predicted curve (black, by this program):  
  
![image](https://github.com/Taomihog/unzipDNA/assets/110962921/710f75ad-8ba1-4234-a182-a5a5bb144cf1)
  
  
For more information on DNA unzipping experiment and theory:  

>[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS.  
>
>[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal.  
>
>[3] Huguet, Josep M., et al. (2010) PNAS.  

