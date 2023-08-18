# DNA Unzipping Curve Calculator
  
## Command For Compilation  

### Compile

Use this command for compiling on both Windows and Linux systems:  
  
>g++ -std=c++20 -fconstexpr-ops-limit=100000000000 main.cpp utils.cpp -DJ_SIZE=200 -DEXT_SIZE=200 -o UnzipDNA.exe

The **`J_SIZE`** macro dictates the lookup table size. Bigger values enhance accuracy but increase compilation times.  
  
For Windows OS, employ the **build_project.bat** script to enable additional options and tests. Successful execution should yield an output like this:  

![image](doc/Compile_time_200x200.png)

*Note: Note: Linux users can translate the **.bat** file to a shell script using ChatGPT*  
  
### Running Tests

Execute the executable, providing the input file name and an optional output file:  
  
>UnzipDNA.exe NEB_H5alpha_Accessory_colonization_factor_AcfD.txt out.csv

The first argument is the input file name. The second argument determines the output file name (optional).  

(DNA sequence is from [genbank CP017100](https://www.ncbi.nlm.nih.gov/nuccore/CP017100))  
  
## Goal of this program

The program calculates unzipping curves rapidly, enabling analysis of numerous genes. Aiming for faster loops, the code precomputes data and stores it.  

I **"constexpr**-ed" most of the calculation overhead. Two look-up tables (LUTs) are created to hold these data. These LUTs are computed and saved in **constexpr std::array** containers (c++20 or above is thus required). The compile time as a drawback, is thousands of times longer than a straightforward C++ program.  

### Multithreading Augmentation

Implemented on Aug/15/2023, multithreading improved execution speed by additional 10-20 times.  

### Performance Milestone

The program now performs over 1,000 times faster than the initial Python code.
  
## DNA unzipping theory

![image](https://github.com/Taomihog/unzipDNA/assets/110962921/710f75ad-8ba1-4234-a182-a5a5bb144cf1)

**Figure above shows DNA unzipping experiment on a 4.4 kb DNA**. The program's predicted unzipping curve (black) aligns well with actual single-molecule measurements (blue), showcasing its accuracy.
  
Further reading on DNA unzipping experiment and theory:  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  
