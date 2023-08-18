#pragma once

#include <string>
#include <vector>

#include "../shared_libs/constants_lib.h"
#include "../shared_libs/math_lib.h"

//a struct to save data point by point
struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

//Calculate DNA sequence's energy
namespace DNAsequence {
    const double EffSaltConc = MyMath::Ln(Condition::SaltConc * 0.001) / 298.0;//298.0 K is where the energy was measured in Huguet paper, it is not equal to Condition::temperature

    const double factor_pNnm = Const::pNnm * Const::Joul / Const::Avogadro / Condition::kT;//convert to 'pN-nm' unit;
    double lookup_bp_energy(char bp1, char bp2);
    std::vector<double> calculate_sequence_energy(const std::string & sequence);
}


std::string creat_path_out(const std::string & path_in);
std::string readtxt_firstline(const std::string & path);


const double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;
dp calculate_array(int extension, const std::vector<double> & seq_energy);