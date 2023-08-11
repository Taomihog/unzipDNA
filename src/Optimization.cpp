#include <functional>
#include <cstdlib>

#include "../include/Optimization.h"

namespace Optimization {
	double Bisection(double x1, double x2, std::function<double(double)> func, double tor) {
		
		if (x1 >= x2) {
			double temp = x1;
			x1 = x2;
			x2 = temp;
		}

		double y1 = func(x1);
		double y2 = func(x2);

		if (y1 * y2 > 0) {
			return inf;//meaning that the root is not found
		}

		double xm, ym;
		int cnt = 0;
		while (cnt < 10000) {
			++cnt;
			if (std::abs(y1) <= tor) {
				return x1;
			}
			if (std::abs(y2) <= tor) {
				return x2;
			}

			xm = (x1 + x2) * 0.5;
			ym = func(xm);

			if (y1 < 0.0 && y2 > 0.0) {
				if (ym > 0) {
					x2 = xm;
					y2 = ym;
				}
				else {
					x1 = xm;
					y1 = ym;
				}
			}
			else if (y1 > 0.0 && y2 < 0.0) {
				if (ym < 0) {
					x2 = xm;
					y2 = ym;
				}
				else {
					x1 = xm;
					y1 = ym;
				}
			}
			else {
				return inf;
			}
		}
		return inf;//meaning that the root is not found
	}
}