#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "../include/DNA_util.h"
#include "../include/Constant_constexpr.h"
#include "../include/LUT_constexpr.h"

namespace DNAsequence {

    void print_LUT_info() {
        std::cout << "LUT extension resolution (nm): " + ext_resolution << std::endl;
        std::cout << "LUT j_unzip resolution: " + j_resolution << std::endl;
    }

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

    double lookup_bp_energy(char bp1, char bp2) {
        int idx1 = Const::bp2idx_map[bp1];
        int idx2 = Const::bp2idx_map[bp2];
        double energy = Const::LUTdH[idx1][idx2] - (Const::LUTm[idx1][idx2] * 
                        Const::EffSaltConc + Const::LUTdS[idx1][idx2] * 0.001 ) * Const::Temperature;
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

namespace LUT_util {

    void test_print_lut (int j0, int k0, int range) {
        int jmax = (j0 + range) > j_size ? j_size : j0 + range;
        int kmax = (k0 + range) > ext_size ? ext_size : k0 + range;

        for (int j = j0; j < jmax; ++j ) {
            for (int k = k0; k < kmax; ++k) {
                std::cout << Lut_force[j][k] << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        for (int j = j0; j < jmax; ++j ) {
            for (int k = k0; k < kmax; ++k) {
                std::cout << Lut_energy[j][k] << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void test_print_lut (int j0, int k0, int range, std::ofstream fout){
        int jmax = (j0 + range) > j_size ? j_size : j0 + range;
        int kmax = (k0 + range) > ext_size ? ext_size : k0 + range;

        for (int j = j0; j < jmax; ++j ) {
            for (int k = k0; k < kmax; ++k) {
                fout << Lut_force[j][k] << '\t';
            }
            fout << std::endl;
        }
    }

    void lookup(double j0, double extension, double & force, double & energy) {
        if (j0 < 0 || j0 >= j_size * j_resolution || extension < 0 || extension >= ext_size * ext_resolution) {
            std::cerr << "Error: lookup_force index out of range";
            return;
        }

        double j = j0 / static_cast<double>(j_resolution);
        double k = extension / static_cast<double>(ext_resolution);

        
        //std::cout << j << '\t' << k << std::endl;
        
        double j1 = std::floor(j);
        double k1 = std::floor(k);

        //test_print_lut(j1, k1, 2);

        double j2 = j1 + 1;
        double k2 = k1 + 1;
        
        //I tried, but the calculation is too tedious...
        //so: https://en.wikipedia.org/wiki/Bilinear_interpolation
        //choose "repeated linear interp"
        //may try polyn fit method later (are they the same?).
        force = Lut_force[j1][k1] * (j2 - j) * (k2 - k) + 
                Lut_force[j1][k2] * (j2 - j) * (k - k1) + 
                Lut_force[j2][k1] * (j - j1) * (k2 - k) + 
                Lut_force[j2][k2] * (j - j1) * (k - k1); 

        energy= Lut_energy[j1][k1] * (j2 - j) * (k2 - k) + 
                Lut_energy[j1][k2] * (j2 - j) * (k - k1) + 
                Lut_energy[j2][k1] * (j - j1) * (k2 - k) + 
                Lut_energy[j2][k2] * (j - j1) * (k - k1); 
        
        //std::cout << force << std::endl << energy << std::endl;
    }
}

namespace {
    //These experiment conditions are hard coded now. backup the original code for future use.
    void read(int argc, char* argv[]) {
        
        if (argc < 2) {
            std::cout << "argv[1]: sequence file location" << std::endl;
            std::cout << "argv[2]: dsDNA arms total length in bp" << std::endl;
            std::cout << "argv[3]: salt concentration [mM]" << std::endl;
            std::cout << "argv[4]: temperature in Kelvin" << std::endl;
            std::cout << "argv[5]: spring constant [pN/nm]" << std::endl;
            return;
        }

        //read the sequence from a file.
        std::string sequence = DNAsequence::readtxt_firstline(argv[1]);
        std::cout << "*Trunk sequence length: \n\t" << sequence.size() << " bp." << std::endl;

        //dsDNA arms' total length.
        int bp_arm = std::stoi(argv[2]);
        std::cout << "*dsDNA arms total length: \n\t" << bp_arm << " bp." << std::endl;

        //get the salt concentration
        double salt = std::stod(argv[3]);
        std::cout << "*Salt concentration: \n\t" << salt << " mM." << std::endl;
        double eff_salt = Math_constexpr::Ln(salt * 0.001) / 298.0;
        std::cout << "*Effective salt concentration: \n\t" << eff_salt << " mM." << std::endl;

        //read the temperature
        double temperature = std::stod(argv[4]);
        std::cout << "*Temperature: \n\t" << temperature << " K." << std::endl;
        double kT = Const::Boltzmann * temperature;

        //read the sprint constant
        double kpillar = std::stod(argv[5]);
        std::cout << "*Spring constant: \n\t" << kpillar << " pN/nm." << std::endl;
        

        //now can calculate the energy to unzip j bp trunck dsDNA, save the data in a vector
        std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
        std::cout << "*Calculated unzipping energy (pNnm): \n\t" ;
        std::for_each(seq_energy.cbegin(), seq_energy.cbegin()+10, [](const double & val) { std::cout << val << ", ";});
        std::cout << "..., " << *(seq_energy.cend()-1) << "\n";

        //look at if the lut is built
        std::cout << "*Test-print the lz_ds: \n\t";
        //code here
        std::cout << std::endl << std::endl;
    }
}
    