#pragma once

#include <vector>
#include <cmath>


namespace {

    constexpr double temperature_measured = 298.0;//the temperature where this set of parameters were measured.

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

    inline int bp2idx(char base) {//todo: need to declare inline explicitly, namespace is not class header!
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

}

inline double calc_eff_salt (double saltConcentration) {
    return std::log(saltConcentration * 0.001) / temperature_measured;
}

inline double lookup_bp_energy(char bp1, char bp2, double eff_salt, double temperature) {
    int idx1 = bp2idx(bp1);
    int idx2 = bp2idx(bp2);
    double energy = LUTdH[idx1][idx2] - (LUTm[idx1][idx2] * eff_salt + LUTdS[idx1][idx2] * 0.001 ) * temperature;

    // cout << LUTdH[idx1][idx2] << '\t' <<
    //         LUTm[idx1][idx2]  << '\t' <<
    //         LUTdS[idx1][idx2]  << endl;

    return energy;
}