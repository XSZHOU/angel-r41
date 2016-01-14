/*
Copyright (c) 2009 Sebastian Steiger, Ratko Veprek, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#include "DomainPoint.h"

using namespace negf;


// =========================================================
// DomainPoint implementations
// =========================================================

/** constructor for 1D domain points */
DomainPoint::DomainPoint(const double& x):
	vertex(-1,x),
	weight(0.0)
{
}

/** constructor for 2D domain points */
DomainPoint::DomainPoint(const double& x, const double& y): 
	vertex(-1,x,y),
	weight(0.0)
{ 
}

/** constructor for 3D domain points */
DomainPoint::DomainPoint(const double& x, const double& y, const double& z): 
	vertex(-1,x,y,z),
	weight(0.0)
{ 
}

double DomainPoint::get_coord_abs() const
{STACK_TRACE(
	double result = 0.0;
	for (uint ii = 0; ii < this->vertex.get_dimension(); ii++) {
		result += vertex.get_coordinate(ii)*vertex.get_coordinate(ii);
	}
	return sqrt(result);
);}


ostream& negf::operator<<(ostream& out, const DomainPoint& point) 
{STACK_TRACE(
	out << "coords = ";
	for (uint ii = 0; ii < point.get_dimension(); ii++) {
		out << point.get_coord(ii) << ", ";
	}
	out <<", weight = " << point.get_weight() << ", index = " << point.get_index();
	return out;  	
);}



// =========================================================
// DomainVertex implementations
// =========================================================

DomainVertex::DomainVertex(int index_, double x_, double y_, double z_):
	dimension(3),
	index(index_)
{
	coord.push_back(x_);
	coord.push_back(y_);
	coord.push_back(z_);
}

DomainVertex::DomainVertex(int index_, double x_, double y_):
	dimension(2),
	index(index_)
{
	coord.push_back(x_);
	coord.push_back(y_);
}

DomainVertex::DomainVertex(int index_, double x_):
	dimension(1),
	index(index_)
{
	coord.push_back(x_);
}

int DomainVertex::get_index() const { 
	return index; 
}

void DomainVertex::set_index(uint index_) { 
	index = index_;
}

/** NEVER change a coordinate during the simulation!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
void DomainVertex::set_coordinate(int xyz, double value)
{STACK_TRACE(
	NEGF_ASSERT(xyz >=0 && xyz < this->dimension, "Invalid dimensionality.");
	this->coord[xyz] = value;
);}

double DomainVertex::get_coordinate(unsigned short int ii) const
{STACK_TRACE(
	NEGF_ASSERT(ii < this->dimension, "Wrong dimensionality.");
	return coord[ii];
);}

double DomainVertex::get_distance_to(DomainVertex * other_vertex) const
{STACK_TRACE(
	NEGF_ASSERT(this->get_dimension()==other_vertex->get_dimension(), "incompatible dimensions.");
	double length2 = 0.0;
	for (uint ii = 0; ii < this->get_dimension(); ii++) {
		double tmp = this->get_coordinate(ii) - other_vertex->get_coordinate(ii);
		length2 += tmp * tmp;
	}
	return std::sqrt(length2);
);}

