# DNA Unzipping Curve Calculator (ver.3)

Build command:

>**g++ -std=c++20 -fconstexpr-ops-limit=100000000000 -o unzip.exe src/main.cpp src/DNA_util.cpp**

*CMAKE somehow only generates build tree for visual studio, but MSVC has issues generating constexpr...*

For more information on DNA unzipping experiment and theory:

>[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS.
>
>[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal.
>
>[3] Huguet, Josep M., et al. (2010) PNAS.

I first wrote the program in Python, and realized that it was too slow for my application. To speed up the calculation, I use look-up table (LUT) to save information that supposes to be calculated in run time. These LUTs are saves in constexpr containers.

It indeed took me some thinking to move majority of calculation overheads from run-time to compile time. After several attempts, I almost **"constexpred"** everyting that can be calculated ahead of time. As a result, this program is 100 times faster than the python program, 6 times faster than the Version 2 program. The only drawback is that the compile time is very long, thousands of times longer than the prototype c++ program (which does calculate many constexpr variables).
