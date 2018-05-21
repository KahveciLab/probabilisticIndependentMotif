/*
 * motif.h
 *
 *  Created on: May 11, 2016
 *      Author: fiona
 */

#ifndef MOTIF_H_
#define MOTIF_H_
#include <vector>
#include <boost/dynamic_bitset.hpp>
using namespace boost; 
using namespace std;

class Motif{
public:
       int id;
	   vector<int> edges;
	   double pri;
	   double priority;
	   double A;
public:
	   Motif(int count): id(count), pri(0), priority(0), A(0){}
	   ~Motif(){
	   	edges.clear();
	   }
};

#endif /* MOTIF_H_ */
