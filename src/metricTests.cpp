#include "../include/metricTests.hpp"

void compareValues(double opti, double brute_pure, double brute_opti,
				   char const*const metric, char const*const dim)
{
	fmt::print("Values for {} {}:\nOptimised:\t{}\nBrute-pure:\t{}\nBrute-opti:\t{}\n", metric, dim, opti, brute_pure,
		   brute_opti);
	if(opti == brute_opti && opti != brute_pure) {
		fmt::print("{}: Difference in psi({}) values, due to ballsize calculation\n", metric, dim);
	}
	if(opti != brute_opti && brute_opti == brute_pure) {
		fmt::print("{}: Difference in psi({}) values, due to count method\n", metric, dim);
	}
	if(opti != brute_opti && brute_opti != brute_pure) {
		fmt::print("{}: Difference in psi({}) values, due to unknown reason\n"
			   "psi({}): {}, psi({}_brute): {}, psi({}_brute-opti): {}\n", metric, dim,
			   dim, opti, dim, brute_pure, dim, brute_opti);
	}
}