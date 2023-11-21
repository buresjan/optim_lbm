// PPolynomial 
// 2015 Radek Fucik
// reads input file defined as
// #1 ... pw polynomial order
// #2 ... number of segments (divisions)
// #3 ... #2+1 rows that describes the divisions
// #4 ... #2x#1 matrix that contains the coefficients of the pw polynomials
#ifndef __PPOLYNOMIAL_H
#define __PPOLYNOMIAL_H

#include "defs.h"

#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include <cstring>
#include <string>

struct PPolynomial
{
	bool m_FLAG_allocated;
	int m_order;  // the order of the polynomials
	int m_pieces; // number of pieces of the polynomials, m_pieces+1 = number of divisions
	double* m_divisions; // array that contains m_pieces+1 increasing points
	double** m_coefs;  // a m_pieces x m_order matrix that contains the coefficients of all polynomials
	int parse_input_file(const char *inputfile);
	int get_index(double x);
	
	bool allocated() { return m_FLAG_allocated; }
	double get(double x);
	double getPeriodic(double x); // get periodic prolong
	
	inline double get_positive(double x) { double ret=get(x); return (ret>0) ? ret : 0; }
	inline double get_negative(double x) { double ret=get(x); return (ret<0) ? ret : 0; }
		
	double get_minx();
	double get_maxx();
	
	PPolynomial(const char *inputfile);
	~PPolynomial();
};

#include "ppolynomial.hpp"

#endif
