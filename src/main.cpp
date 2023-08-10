#include <vector>
#include <array>
#include <iostream>
#include <algorithm>

#include "DNAmech.h"
#include "DNAbp.h"

using namespace std;

static_assert(lz_ds[lsize-1]);

int main() {
    for_each(lz_ds.cbegin(), lz_ds.cend(), [](const double & val) { cout << val << '\t';});
    cout << endl << endl;

    cout << DNAbp::BPEnergy(1, 2, 0.01, 298.0) << endl;
    cout << endl << endl;
    return 0;

}