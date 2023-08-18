#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

#include "include/utils.h"
#include "include/lut_lib.h"

//=============================================util functions==========================================

//single read. todo: read the fasta directly, multiple line read.
std::string readtxt_firstline(const std::string & path) {
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

namespace DNAsequence {

    double lookup_bp_energy(char bp1, char bp2) {
        int idx1 = BPEnergy::bp2idx_map[bp1];
        int idx2 = BPEnergy::bp2idx_map[bp2];
        double energy = BPEnergy::LUTdH[idx1][idx2] - (BPEnergy::LUTm[idx1][idx2] * 
                        EffSaltConc + BPEnergy::LUTdS[idx1][idx2] * 0.001 ) * Condition::Temperature;
        return - energy * factor_pNnm;//convert to 'pN-nm' unit;
    }

    std::vector<double> calculate_sequence_energy(const std::string & sequence) {
        //from the DNA sequence, calculate the energy at every j_unzipped.

        std::vector<double> sequenceEnergy;
        if (sequence.size() < 2) {
            std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
            return sequenceEnergy;
        }

        double accum = 0.0;
        *std::back_inserter(sequenceEnergy) = 0.0;
        std::transform(sequence.cbegin() + 1, sequence.cend(), std::back_inserter(sequenceEnergy), 
        [&accum](const char& bp2) {
            accum += lookup_bp_energy(*(&bp2 - 1), bp2);
            return accum;
            }
        );//equivalent to std::partial_sum()

        return sequenceEnergy;
    }

}

dp calculate_array(int extension, const std::vector<double> & seq_energy) {
        
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        temp_f.at(j) = Force(j,extension);
        temp_e.at(j) = Energy(j,extension);
        // temp_f.at(j) = lookup(j,extension,Lut_force);
        // temp_e.at(j) = lookup(j,extension,Lut_energy);
        
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
