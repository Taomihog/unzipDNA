#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>

//what needs to be included at the end
#include "../include/LUT_constexpr.h"
#include "../include/DNA_util.h"
//Optimization.h is removed from this project because <functional> header cannot be used to in constexpr functions

//=============================================code==========================================

const double energy_threshold = 100.0;//don't calculate exp(-e/kT) and set probability to 0;
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
    std::ofstream fout(path_out + "out_" + argv[1]);
    if (!fout) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }
    fout << "extension (nm),average force (pN),sd force (pN),average bp unzipped, sd bp unzipped" << std::endl;

    //allocate some memory
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    //loop start
    for (int extension = 0; extension < 1.2 * seq_energy.size() + 0.35 * Const::ArmLength; ++extension) {

        double force, energy;

        double min_e = 1.0e100;

        for (int j = 0; j < seq_energy.size(); ++j) {
            LUT_util::lookup(j,extension,force,energy);
            temp_f.at(j) = force;
            temp_e.at(j) = energy + seq_energy[j];
            if (min_e > temp_e.at(j)) {
                min_e = temp_e.at(j);
            }
            //std::cout << j << '\t' << temp_f.at(j) << '\t' << temp_e.at(j) << std::endl;
        }

        std::for_each(temp_e.begin(), temp_e.end(), [min_e](double &val) {
            val -= min_e;
            if (val > energy_threshold) {
                val = 0.0;
            } else {
                val = std::exp(-val/Const::kT);
            }
        });
        
        double prob = 0;
        double Fprob = 0;
        double FFprob = 0;
        double Jprob = 0;
        double JJprob = 0;
        for (int j = 0; j < seq_energy.size(); ++j) {
            prob += temp_e.at(j);
            Fprob += temp_f.at(j) * temp_e.at(j);
            FFprob += temp_f.at(j) * temp_f.at(j) * temp_e.at(j);
            Jprob += j * temp_e.at(j);
            JJprob += j * j * temp_e.at(j);
        }
        
        //output
        char delimit = '\t';
        if (extension % 1000 == 0) {
            std::cout << extension << delimit << Fprob/prob << delimit << std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob)) 
                                    << delimit << Jprob/prob << delimit << std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob)) << std::endl;
        }
        delimit = ',';
        fout << extension << delimit << Fprob/prob << delimit << std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob)) 
                                << delimit << Jprob/prob << delimit << std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob)) << std::endl;
    }

    // Close the file
    fout.close();
    std::cout << "Table has been written to 'table.txt'." << std::endl;
    
    return 0;
}