#include <iostream>
#include <fstream>

#include "LUT_constexpr.h"

namespace UnzipLUT {

    void print_LUT_info() {
        std::cout << "LUT extension resolution (nm): " + ext_resolution << std::endl;
        std::cout << "LUT j_unzip resolution: " + j_resolution << std::endl;
    }

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
            force = - 1.0;
            energy = -1.0;
            return;
        }

        double j = j0 / static_cast<double>(j_resolution);
        double k = extension / static_cast<double>(ext_resolution);
        
        //because j and k are positive.
        //static_cast<int> has the same result as std::floor() from cmath lib
        //Legend says that cast is 3 times faster
        double j1 = static_cast<int>(j);
        double k1 = static_cast<int>(k);

        //test_print_lut(j1, k1, 2);

        double j2 = j1 + 1;
        double k2 = k1 + 1;

        //https://en.wikipedia.org/wiki/Bilinear_interpolation

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