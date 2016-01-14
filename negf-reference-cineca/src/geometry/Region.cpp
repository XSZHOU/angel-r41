/*
Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#include "Region.h"
#include <stdio.h>
#include <vector>


using namespace negf;

	
Region::Region(const char* name_): 
	name(name_),
	material_name("unknown_material"),
	material_molefraction(-1.0),
	num_dfise_vertices(0)
{
	this->mat = 0;	
	this->elements.clear();
}
Region::Region(const string& name_):
	name(name_),
	material_name("unknown_material"),
	material_molefraction(-1.0),
	num_dfise_vertices(0)
{
	this->mat = 0;
	this->elements.clear();
}

Region::~Region() {
	this->mat = 0;	
}

const string& Region::get_name() const {
	return this->name;	
}


/*const*/PropertyContainer<double>* Region::get_material() const 
{STACK_TRACE(
	NEGF_FASSERT(this->mat != 0, "no material was assigned to region %s", this->get_name().c_str());
	return this->mat;	
);}

const string& Region::get_material_name() const {
	return this->material_name;
}

double Region::get_material_molefraction() const 
{STACK_TRACE(
	NEGF_ASSERT(material_molefraction!=-1.0, "molefraction was not yet assigned!");
	return this->material_molefraction;
);}

bool Region::has_molefraction() const 
{STACK_TRACE(
	if (fabs(this->material_molefraction + 1.0) < 1e-10)
		return false;
	else
		return true;
);}

void Region::set_material(PropertyContainer<double>* mat_) 
{STACK_TRACE(
	NEGF_ASSERT(mat_!=NULL, "null pointer encountered.");
	if (this->material_name!="unknown_material" && this->material_name!=mat_->get_name()) {
		logmsg->emit(LOG_WARN,"Replacing material \"%s\" of region \"%s\" with \"%s\"",
			this->material_name.c_str(), this->get_name().c_str(), mat_->get_name().c_str());
		// do not throw error because this is regularly the case with ternary materials!
		// AlGaAs --> AlGaAs0.3
	}
	this->mat = mat_;
	this->material_name = mat->get_name();
	// molefraction is NOT assigned!
);}

void Region::set_material_name(const char* key) 
{STACK_TRACE(
	this->material_name = key;
	// do not display warning message because at the time the geometry is creted the name is everything that is known of the material!
	/*if (this->mat==NULL) {
		logmsg->emit(LOG_WARN,"warning: a material name was assigned to region \"%s\" but not the material itself.",
				this->get_name().c_str());
	}*/
);}

void Region::set_material_molefraction(const double molefraction) 
{STACK_TRACE(
	NEGF_ASSERT(material_molefraction==-1.0, "molefraction was already assigned!");
	this->material_molefraction = molefraction;	
);}


uint Region::get_index() const
{STACK_TRACE(
	if (index_ready) 
		return index; 
	else 
		NEGF_EXCEPTION("Index for this region has not been set yet."); 
	return 0; // never reached
);}


void Region::add_element(Element * elem)
{STACK_TRACE(
	this->elements.push_back(elem);
);}


void Region::set_dfise_region_vertex_numbers(const vector<uint> & vertex_indices)
{STACK_TRACE(
	this->num_dfise_vertices = vertex_indices.size();
	this->dfise_vertices = vertex_indices;
);}

bool Region::is_oxide() const 
{STACK_TRACE(
	bool oxide = false;
	for (uint nn=0; nn < Constants.oxide_names.size(); nn++) {
		if (this->get_material_name()==Constants.oxide_names[nn]) {	
			oxide = true; 
			break;	
		}
	}
	return oxide;
);}
