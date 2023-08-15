#include <iostream>
#include <fstream>

namespace UnzipLUT {
    void print_LUT_info();
    void test_print_lut (int j0, int k0, int range);
    void test_print_lut (int j0, int k0, int range, std::ofstream fout);
    void lookup(double j0, double extension, double & force, double & energy);
}
