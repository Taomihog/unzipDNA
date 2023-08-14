#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "../include/LUT_constexpr.h"
#include "../include/DNA_util.h"
//Optimization.h has been removed from this project, because <functional> header cannot be used in constexpr functions

//=============================================code==========================================

const double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;
int main(int argc, char * argv[]) {

    if (argc <2) {
        std::cerr << "Please provide a filename (argv[1])";
        return 1;
    }
    std::string sequence = DNAsequence::readtxt_firstline(argv[1]);
    std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
    //test, the "41740.760955375" is from my python code
    std::cout << "Sequence energy calculate difference from python code: "<< *(seq_energy.cend() - 1) - 41740.760955375 << std::endl;// no difference, good
    
    //prepare file for output
    std::string path_out;
    path_out = path_out + "out_" + argv[1]; 
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
    for (int extension = 0; extension < 1.2 * seq_energy.size() + 0.35 * Const::ArmLength; ++extension) {

        double min_e = 1.0e100;
        for (int j = 0; j < seq_energy.size(); ++j) {
            
            LUT_util::lookup(j,extension,temp_f.at(j),temp_e.at(j));
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

            //slightly faster than code below:
            // temp_e.at(j) = temp_e.at(j) - min_e;
            // temp_e.at(j) = temp_e.at(j) > energy_threshold ? 0.0 : std::exp(-temp_e.at(j));
            
            // prob += temp_e.at(j);
            // Fprob += temp_f.at(j) * temp_e.at(j);
            //FFprob += temp_f.at(j) * temp_f.at(j) * temp_e.at(j);
            // Jprob += j * temp_e.at(j);
            // JJprob += j * j * temp_e.at(j);
        }
        
        //output
        double dna_length = extension - Fprob/prob/Const::PillarStiffness;//delete pillar deformation from the extension

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