#include "../include/defs.hpp"

double angleDist(double i, double j) {
	return dwrap(std::fabs(i - j), TWOPI);
}

double regularDist(double i, double j) {
	return std::fabs(i - j);
}