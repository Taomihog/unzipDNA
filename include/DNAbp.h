#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <limits>

using namespace std;

//file read to be implemented
std::string defaultSequence = { "CGCAATCATCGATATCTCGTAATCACGTGCAAGGCCTACTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCG" };//example data

namespace {
    constexpr double temperature_measured = 298;//the temperature where this set of parameters were measured.
}

namespace DNAbp {
    constexpr int LUTsize = 4;
    constexpr double LUTdH[LUTsize][LUTsize] = {
        {-7.28, -4.63, -5.21, -5.80},//aa, at, ag, ac
        {-8.31, -7.28, -8.96, -8.16},//ta, tt, tg, tc
        {-8.16, -5.80, -8.57, -10.1},//ga, gt, gg, gc
        {-8.96, -5.21, -9.66, -8.57} //ca, ct, cg, cc
    };
    constexpr double LUTdS[LUTsize][LUTsize] = {
        {-20.28, -11.62, -12.89, -14.46},//aa, at, ag, ac
        {-25.06, -20.28, -24.48, -22.46},//ta, tt, tg, tc
        {-22.46, -14.46, -22.30, -25.96},//ga, gt, gg, gc
        {-24.48, -12.89, -24.43, -22.30} //ca, ct, cg, cc
    };
    constexpr double LUTm[LUTsize][LUTsize] = {
        {0.145, 0.117, 0.070, 0.099},//aa, at, ag, ac
        {0.091, 0.145, 0.091, 0.155},//ta, tt, tg, tc
        {0.155, 0.099, 0.063, 0.079},//ga, gt, gg, gc
        {0.091, 0.070, 0.132, 0.063} //ca, ct, cg, cc
    };

    inline int bp2idx(char base) {//todo: need to declare inline otherwise will be an error somehow!
        switch (base) {
        case 'a':
        case 'A':
            return 0;
        case 't':
        case 'T':
            return 1;
        case 'g':
        case 'G':
            return 2;
        case 'c':
        case 'C':
            return 3;
        default:
            return -1;
        }
    };

    inline double BPEnergy(int idx1, int idx2, double eff_salt, double temperature) {
        return LUTdH[idx1][idx2] - 
            LUTm[idx1][idx2] * eff_salt * temperature / temperature_measured +
            LUTdS[idx1][idx2] * 0.001 * temperature;
    }
}


//data in original form:
//bpLUT = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT" };
//LUTdH = { -7.28, -5.8, -5.21, -4.63, -8.96, -8.57, -9.66, -5.21, -8.16, -10.1, -8.57, -5.8, -8.31, -8.16, -8.96, -7.28 };
//LUTdS = { -20.28, -14.46, -12.89, -11.62, -24.48, -22.3, -24.43, -12.89, -22.46, -25.96, -22.3, -14.46, -25.06, -22.46, -24.48, -20.28 };
//LUTm  = { 0.145, 0.099, 0.07, 0.117, 0.091, 0.063, 0.132, 0.07, 0.155, 0.079, 0.063, 0.099, 0.091, 0.155, 0.091, 0.145 };