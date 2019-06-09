#ifndef NN_HPP
#define NN_HPP

#include "../include/vpUtil.hpp"
#include "../include/vpTree.hpp"
#include <vector>

/**
 * Finds K nearest neighbours to each point i in points.
 * points is a 2 dimensional point cloud in domain [0, L]
 * Wrapping indicates whether neighbour proximitiy
 * wraps around border.
 * neighbours is output 1D array. K nearest neighbours for
 * first particle are stored in first K elements, neighbours
 * for next particle are in next K elements, etc.
 * That is, the kth nearest neighbour to particle i
 * will be stored in [i*K+k]
 */
void getNearestNeighbours(std::vector<VPoint<2> > &points,
						  int const N,
						  int const K,
						  double const L,
						  int * neighbours);

/**
* \brief Finds the k nearest points
* This function is hard-coded in the linear space.
* Is only used for finding neighbours for the GTE algorithm, to determine which particles to send to gomez.
* Not used in any gomez functions
* \param dataPts - Point cloud to search
* \param point - Query point
* \param numPts - Number of points in point cloud
* \param k - Number of neighbours to search for
* \param dim - Dimensions in point cloud
* \param linearSize - Wrapping parameter
* \return List of nearest neighbours
*/
std::vector<int> brutedistkPoints(std::vector<VPoint<2> > const & dataPts, int const point,
								  int const numPts, int const k, int const dim, 
								  double linearSize);


#endif //NN_HPP