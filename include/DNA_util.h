#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include "Constant_constexpr.h"

//DNA sequence related stuff, include DNA sequence io, 
//DNA sequence's energy,
//read an parse DNA sequence data,//etc.

namespace DNAsequence {
    //single read. todo: read the fasta directly, multiple line read.
    std::string readtxt_firstline(const std::string & path);

    //utilities to convert bp to energy
    const double factor_pNnm = Const::pNnm * Const::Joul / Const::Avogadro;//convert to 'pN-nm' unit;
    double lookup_bp_energy(char bp1, char bp2);
    std::vector<double> calculate_sequence_energy(const std::string & sequence);
}

namespace LUT_util {
    void print_LUT_info();
    void test_print_lut (int j0, int k0, int range);
    void test_print_lut (int j0, int k0, int range, std::ofstream fout);
    void lookup(double j0, double extension, double & force, double & energy);
}
