#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "../include/DNAmech.h" //LUTs: lforce, lz_ss, lz_ds, le_ss, le_ds
#include "../include/DNAbp.h"
#include "../include/Constant.h"
#include "../include/Optimization.h"

//#define test

std::string readtxt_firstline(const std::string & path) {
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return "0";
    }

    std::string line;
    if (getline(file, line)) {
        std::cout << "Trunk sequence read successfully." << std::endl;
    } else {
        std::cerr << "Failed to read a line from the file." << std::endl;
    }

    file.close();

    return line;
}

std::vector<double> calculate_sequence_energy(const std::string & sequence, double eff_salt, double temperature) {
    //from the DNA sequence, calculate the energy at every j_unzipped.

    std::vector<double> sequenceEnergy;
    if (sequence.size() < 2) {
        std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
        return sequenceEnergy;
    }

    double factor = - Constant::pNnm * Constant::Joul / Constant::Avogadro;//convert to 'pN-nm' unit

    //equivalent to std::partial_sum()
    double accum = 0.0;
    *std::back_inserter(sequenceEnergy) = 0.0;
    std::transform(sequence.cbegin() + 1, sequence.cend(), std::back_inserter(sequenceEnergy), 
    [factor, eff_salt, temperature, &accum](const char& bp2) {
        accum += factor * lookup_bp_energy(*(&bp2 - 1), bp2, eff_salt, temperature);
        return accum;
        }
    );

    return sequenceEnergy;
}

double calc_force (double extension, double bp_arm, double kpillar, double kT, const std::vector<double> & seq_energy) {

    //define range of partition function calculation
    const double force_upper = 25.0;
    const double force_lower = 5.0;
    
    //calculate the jmax and jmin at force_max and force_min\
    //be carefull that when j decreases when force increases
    #ifdef test
    int jmax =     std::floor((extension - force_lower/kpillar - bp_arm * lz_ds_f(force_lower)) / lz_ss_f(force_lower)) / 2.0;
    int jmin = 1 + std::floor((extension - force_upper/kpillar - bp_arm * lz_ds_f(force_upper)) / lz_ss_f(force_upper)) / 2.0;
    #else
    int jmax =     std::floor((extension - force_lower/kpillar - bp_arm * lookup(lz_ds, force_lower)) / lookup(lz_ss,force_lower)) / 2.0;
    int jmin = 1 + std::floor((extension - force_upper/kpillar - bp_arm * lookup(lz_ds, force_upper)) / lookup(lz_ss,force_upper)) / 2.0;
    #endif
    //coerce jmax and jmin
    jmax = jmax > seq_energy.size() ? seq_energy.size() -1 : jmax;
    jmin = jmin < 0 ? 0 : jmin;
    //std::cout << "jmin: " << jmin << ", jmax: " << jmax << "\t";

    //calculate the energy, force, etc at each j between jmin and jmax.
    std::vector<double> force_arr;
    std::vector<double> energy_arr;
    std::vector<double> j_arr;
    double energy_min = 1.0e10;
    for (int j = jmin; j <= jmax; ++j) {

        j_arr.emplace_back(j);
        //std::cout << "emplace_back: j: " << j;

        //binary search to get the force 
        double force = Optimization::Bisection(force_min, force_max, [=](double f){
            #ifdef test
            return extension - f / kpillar - bp_arm * lz_ds_f(f) - 2.0 * j * lz_ss_f(f);
            #else
            return extension - f / kpillar - bp_arm * lookup(lz_ds, f) - 2.0 * j * lookup(lz_ss,f);
            #endif
        });

        force_arr.emplace_back(force);//which is faster? back_inserter vs emplace_back
        //std::cout << ", force: " << force;

        //calculate the total energy
        #ifdef test
        double energy = 0.5 * force * force / kpillar + bp_arm * le_ds_f(force) + 
                        2.0 * j * le_ss_f(force) + seq_energy[j];
        #else
        double energy = 0.5 * force * force / kpillar + bp_arm * lookup(le_ds,force) + 
                        2.0 * j * lookup(le_ss, force) + seq_energy[j];
        #endif

        energy_arr.emplace_back(energy);
        //std::cout << ", energy: " << energy << ".\n";

        //record the min energy
        if (energy < energy_min) {
            energy_min = energy;
        }
    }
    std::vector<double> probability_arr;
    std::transform(energy_arr.cbegin(), energy_arr.cend(), std::back_inserter(probability_arr), [energy_min, kT](const double & energy) {
        return exp((energy_min - energy) / kT);
    });

    if (probability_arr.size() != force_arr.size()) {
        std::cerr << "Error: somehow arrays sizes are not equal.";
    }

    double p_sum = 0.0;
    double pf_sum = 0.0;
    double pj_sum = 0.0; 
    for (int i = 0; i < probability_arr.size(); ++i) {
        p_sum += probability_arr.at(i);
        pf_sum += probability_arr.at(i) * force_arr.at(i);
        pj_sum += probability_arr.at(i) * j_arr.at(i);
    }

    double f_avg = pf_sum / p_sum;

    return f_avg;
}

//command to build and run
//g++ -std=c++20 -o a.exe src/main.cpp src/Optimization.cpp && a.exe exampleData.txt 2200 100 298.15 0.07406
int main(int argc, char * argv[]) {

    //============================================read parameters======================================
    if (argc < 6) {
        std::cout << "argv[1]: sequence file location" << std::endl;
        std::cout << "argv[2]: dsDNA arms total length in bp" << std::endl;
        std::cout << "argv[3]: salt concentration [mM]" << std::endl;
        std::cout << "argv[4]: temperature in Kelvin" << std::endl;
        std::cout << "argv[5]: spring constant [pN/nm]" << std::endl;
        return 0;
    }

    //read the sequence from a file.
    std::string sequence = readtxt_firstline(argv[1]);
    std::cout << "*Trunk sequence length: \n\t" << sequence.size() << " bp." << std::endl;

    //dsDNA arms' total length.
    int bp_arm = std::stoi(argv[2]);
    std::cout << "*dsDNA arms total length: \n\t" << bp_arm << " bp." << std::endl;

    //get the salt concentration
    double salt = std::stod(argv[3]);
    std::cout << "*Salt concentration: \n\t" << salt << " mM." << std::endl;
    double eff_salt = calc_eff_salt(salt);
    std::cout << "*Effective salt concentration: \n\t" << eff_salt << " mM." << std::endl;

    //read the temperature
    double temperature = std::stod(argv[4]);
    std::cout << "*Temperature: \n\t" << temperature << " K." << std::endl;
    double kT = Constant::Boltzmann * temperature;

    //read the sprint constant
    double kpillar = std::stod(argv[5]);
    std::cout << "*Spring constant: \n\t" << kpillar << " pN/nm." << std::endl;

    //now can calculate the energy to unzip j bp trunck dsDNA, save the data in a vector
    std::vector<double> seq_energy = calculate_sequence_energy(sequence, eff_salt, temperature);
    std::cout << "*Calculated unzipping energy (pNnm): \n\t" ;
    std::for_each(seq_energy.cbegin(), seq_energy.cbegin()+10, [](const double & val) { std::cout << val << ", ";});
    std::cout << "..., " << *(seq_energy.cend()-1) << "\n";

    //look at if the lut is built
    std::cout << "*Test-print the lz_ds: \n\t";
    std::for_each(lz_ds.cbegin(), lz_ds.cbegin()+7, [](const double & val) { std::cout << val << ", ";});
    std::cout << std::endl << std::endl;

    //========================================calculate unzipping curve===================================
    for (double extension = 0.0; extension < 6000.0; extension += 1.0) {
        double force = calc_force (extension, bp_arm, kpillar, kT, seq_energy);
        std::cout << "force is " << force << " pN at extension = " << extension- force/kpillar << " nm. \n";
    }
    std::cout << std::endl << std::endl;

    
    return 0;

}