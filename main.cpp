#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <chrono>

#include "shared_libs/ThreadPool.h"//https://github.com/progschj/ThreadPool
#include "include/lut_lib.h"//LUT library
#include "include/utils.h"//DNA unzip energy calculation library

//==========================================what is DNA unzipping experiment==============================================
//the DNA structure used in an unzipping experiment is something like this:
//
//            Force
//             ↑
//             ||
//             ||
//               =================
//             ||
//             ||
//             ↓
//            Force
//
//This is called a "Y structure". A Y structure has 2 dsDNA arms and a dsDNA trunk.
//when the two arms are stretched as shown, the trunk will be unzipped (ie, paired nucleotides will be separated).
//the force during this "unzip" precess is measured against the end-to-end distance of the arms. 
//The sequence of the trunk determines the "force vs extension" profile.
//AND can be calculated theoretically.

//for more information, see:

//[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS.
//[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal.
//[3] Huguet, Josep M., et al. (2010) PNAS.


// Define an alias for the high-resolution clock
using Clock = std::chrono::high_resolution_clock;

//==================================================main===========================================

int main(int argc, char * argv[]) {

    auto start = Clock::now();

    if (argc <2) {
        std::cerr << "Please provide a filename (argv[1])";
        return 1;
    }
    
    const std::string sequence = readtxt_firstline(argv[1]);
    const std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
    //test, the "41740.760955375" is from my python code
    //std::cout << "Sequence energy differs from python program by "<< *(seq_energy.cend() - 1) - 10140.0933068047 << std::endl;// no difference, good
    
    int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads: " << numThreads << std::endl;
    ThreadPool pool(numThreads);

    std::vector< std::future< dp > > results;

    for(int extension = 1; extension < static_cast<int>(1.2 * seq_energy.size()); ++extension) {
        results.emplace_back(
            pool.enqueue([extension, &seq_energy]{
                return calculate_array(extension, seq_energy);
            })
        );
    }

    std::vector<dp> result_array;
    for(auto && result: results)
        result_array.emplace_back(result.get());

    //prepare file for output
    const std::string path_out = argc >= 3 ? argv[2] : creat_path_out(argv[1]);

    std::ofstream fout(path_out);
    if (!fout) {
        std::cerr << "Error opening file. Cannot save result" << std::endl;
        return 1;
    }

    fout << "total extension (nm),DNA extension (nm),average force (pN),sd force (pN),average bp unzipped,sd bp unzipped" << std::endl;

    char delimit = ',';
    std::for_each(result_array.cbegin(), result_array.cend(), [&fout, delimit](const dp & point) {
        fout << point.extension_total << delimit;
        fout << point.extension_DNA << delimit;
        fout << point.force_average << delimit;
        fout << point.force_SD << delimit;
        fout << point.junzipped_average << delimit;
        fout << point.junzipped_SD << std::endl;;
    });

    // Close the file
    fout.close();
    std::cout << "Result has been written to '" + path_out << "'." <<std::endl;
    
    auto end = Clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;
    
    return 0;
}