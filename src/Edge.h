/*
 * edge.h
 *
 *  Created on: May 11, 2016
 *      Author: fiona
 */

#ifndef EDGE_H_
#define EDGE_H_
#include <armadillo>
using namespace arma;

class Edge{
public:
       mat temp; // record the index that can complete a motif and also its template that can be completed
       colvec degree; // record the index that indicates a motif
};

#endif /* Edge_H_ */
