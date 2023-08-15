#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

//need the condition from Constant_constexpr.h
#include "../include/Constant_constexpr.h"

//Calculate DNA sequence's energy
namespace DNAsequence {
    const double factor_pNnm = Const::pNnm * Const::Joul / Const::Avogadro / Condition::kT;//convert to 'pN-nm' unit;
    double lookup_bp_energy(char bp1, char bp2);
    std::vector<double> calculate_sequence_energy(const std::string & sequence);
}