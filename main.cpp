#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

//main function still need values in this header
#include "include/Constant_constexpr.h"
//LUT library
#include "UnzipLUT/UnzipLUT.h"
//DNA unzip energy calculation library
#include "BPEnergy/BPEnergy.h"
//Optimization.h has been removed from this project, because <functional> header cannot be used in constexpr functions

constexpr double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;

//=============================================code==========================================


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


int main(int argc, char * argv[]) {

    if (argc <2) {
        std::cerr << "Please provide a filename (argv[1])";
        return 1;
    }
    std::string sequence = readtxt_firstline(argv[1]);
    std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
    //test, the "41740.760955375" is from my python code
    std::cout << "Sequence energy calculate difference from python code: "<< *(seq_energy.cend() - 1) - 10140.0933068047 << std::endl;// no difference, good
    
    //prepare file for output
    std::string path_in(argv[1]);

    size_t lastSlashIndex = path_in.rfind('/');
    std::cout << lastSlashIndex << std::endl;
    std::string parentPath;
    if (lastSlashIndex != std::string::npos) {
        std::cout << "last slash index is not zero" <<std::endl;
        parentPath = path_in.substr(0, lastSlashIndex) + "/";
    }
    std::string fullname = path_in.substr(lastSlashIndex + 1, std::string::npos);
    
    std::string fileExtension;
    size_t dotIndex = fullname.rfind('.');
    std::cout << "parent path: " << parentPath << std::endl;
    std::string path_out = parentPath + "out_" + fullname.substr(0, dotIndex) + ".csv";

    std::cout << "path out: " << path_out << std::endl;

    std::ofstream fout(path_out);
    if (!fout) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }
    fout << "extension (nm),average force (pN),sd force (pN),average bp unzipped,sd bp unzipped" << std::endl;

    std::cout << "extension (nm)\taverage force (pN)\tsd force (pN)\taverage bp unzipped\tsd bp unzipped" << std::endl;

    //allocate some memory for calculation
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    //loop start
    for (int extension = 0; extension < 1.2 * seq_energy.size() + 0.35 * Condition::ArmLength; ++extension) {

        double min_e = 1.0e100;
        for (int j = 0; j < seq_energy.size(); ++j) {
            
            UnzipLUT::lookup(j,extension,temp_f.at(j),temp_e.at(j));
            temp_e.at(j) +=seq_energy[j];

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
        
        //output
        double dna_length = extension - Fprob/prob/Condition::PillarStiffness;//delete pillar deformation from the extension

        char delimit = '\t';
        if (static_cast<int>(extension) % 1000 == 0) {
            std::cout << dna_length << delimit << Fprob/prob << delimit << std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob)) 
                                    << delimit << Jprob/prob << delimit << std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob)) << std::endl;
        }
        delimit = ',';
        fout << dna_length << delimit << Fprob/prob << delimit << std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob)) 
                                << delimit << Jprob/prob << delimit << std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob)) << std::endl;

    }

    // Close the file
    fout.close();
    std::cout << "Result has been written to '" + path_out << "'." <<std::endl;
    
    return 0;
}