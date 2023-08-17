#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <chrono>

#include "include/ThreadPool.h"//https://github.com/progschj/ThreadPool
#include "include/Constant_constexpr.h"//main function still need values in this header
#include "UnzipLUT/UnzipLUT.h"//LUT library
#include "BPEnergy/BPEnergy.h"//DNA unzip energy calculation library
//Optimization.h has been removed from this project, because <functional> header cannot be used in constexpr functions, has to hard-code the function



// Define an alias for the high-resolution clock
using Clock = std::chrono::high_resolution_clock;

constexpr double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;


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



//=============================================util functions==========================================

//single read. todo: read the fasta directly, multiple line read.
std::string readtxt_firstline(const std::string & path){
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return "0";
    }

    std::string line;
    if (getline(file, line)) {
        std::cout << "Trunk sequence '" + path + "' read successfully." << std::endl;
    } else {
        std::cerr << "Failed to read a line from the file." << std::endl;
    }

    file.close();

    return line;
}

//generate a file path for output result from the input file path.
std::string creat_path_out(const std::string & path_in)
{
    size_t lastSlashIndex = path_in.rfind('/');

    std::string parentPath;
    if (lastSlashIndex != std::string::npos) {
        parentPath = path_in.substr(0, lastSlashIndex) + "/";
    }

    std::string fullname = path_in.substr(lastSlashIndex + 1, std::string::npos);

    return parentPath + "out_" + fullname.substr(0, fullname.rfind('.')) + ".csv";;
}

//a struct to save data point by point
struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

dp calculate_array(int extension, const std::vector<double> & seq_energy) {
        
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        
        UnzipLUT::lookup(j,extension,temp_f.at(j),temp_e.at(j));
        temp_e.at(j) += seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    double prob = 0;
    double Fprob = 0;
    double FFprob = 0;
    double Jprob = 0;
    double JJprob = 0;
    double temp1,temp2,temp3;
    for (int j = 0; j < seq_energy.size(); ++j) {

        temp1 = temp_e.at(j) - min_e;
        temp1 = temp1 > energy_threshold ? 0.0 : std::exp(-temp1);
        temp2 = temp_f.at(j);

        prob += temp1;

        temp3 = temp2 * temp1;
        Fprob += temp3;
        FFprob += temp3 * temp2;

        temp3 = j * temp1;
        Jprob += temp3;
        JJprob += j * temp3;

    }

    dp point;
    point.extension_total = extension;
    point.force_average = Fprob/prob;
    point.extension_DNA = extension - point.force_average/Condition::PillarStiffness;
    point.force_SD = std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob));
    point.junzipped_average = Jprob/prob;
    point.junzipped_SD = std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob));
    return point;
}



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
    
    //const std::string path_in(argv[1]);
    std::string path_out = argc >= 3 ? argv[2] : creat_path_out(argv[1]);

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