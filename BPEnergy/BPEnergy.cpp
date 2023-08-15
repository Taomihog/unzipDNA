#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

#include "BPEnergy.h"
#include "BPEnergy_constexpr.h"

namespace DNAsequence {

    double lookup_bp_energy(char bp1, char bp2) {
        int idx1 = BPEnergy::bp2idx_map[bp1];
        int idx2 = BPEnergy::bp2idx_map[bp2];
        double energy = BPEnergy::LUTdH[idx1][idx2] - (BPEnergy::LUTm[idx1][idx2] * 
                        BPEnergy::EffSaltConc + BPEnergy::LUTdS[idx1][idx2] * 0.001 ) * Condition::Temperature;
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