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
#ifndef MATERIALDATABASE_H_NEGF
#define MATERIALDATABASE_H_NEGF

#include "all.h"

#include "PropertyContainer.h"
#include "TernaryPropertyContainer.h"

using namespace std;

namespace negf
{		
	/** material properties container
	* 
	* the material properties container holds the information on all materials avaiable.
	* materials in the calculation are defined by their name (identification string). if 
	* a certain property of a material referenced in the grid file is requested, the 
	* material database tries to find it internally. if it is not yet loaded, it
	* performs a lookup in all path directories to locate a file called
	* <identificationstring>.mat which should hold the material properties. 
	*/
	class MaterialDatabase
	{
	public:
	
		typedef map<string,	PropertyContainer<double>*>::iterator material_iterator;
		typedef map<string,	PropertyContainer<double>*>::const_iterator material_citerator;
		typedef map<string,	TernaryPropertyContainer<double>*>::iterator ternary_material_iterator;
		typedef map<string,	TernaryPropertyContainer<double>*>::const_iterator ternary_material_citerator;
	
		MaterialDatabase() throw (Exception *);
		~MaterialDatabase() throw (Exception *);
		
		void   add_search_path		(const char* path) 					throw (Exception *);
		double get					(const char* mat, const char* key) 	throw (Exception *);
		bool   material_exists		(const char* key) const				throw (Exception *);
		bool   ternary_material_exists(const char* key) const 			throw (Exception *);
		bool   material_is_valid	(const char* key) 					throw (Exception *);
		bool   property_is_set		(const char* mat, const char* key) 	throw (Exception *);
		void   load_material		(const char* name) 					throw (Exception *);	
		void   load_ternary_material(const char* name, const double molefraction) throw (Exception *);	
		void   add_material			(const string& name, PropertyContainer<double>* material) throw (Exception *);
		void   add_ternary_material(const string& name, TernaryPropertyContainer<double>* prop)  throw (Exception *);
		
		PropertyContainer<double>* get_material(const char* key) const throw (Exception *); 
		PropertyContainer<double>* get_material(const int &id)   const throw (Exception *);	
		TernaryPropertyContainer<double>* get_ternary_material(const char* key) const throw (Exception *); 
		size_t get_num_materials() const { return this->materials.size(); }
		
								
	protected:
		list<string> paths;
		map<string,	PropertyContainer<double>*>			materials;
		map<string,	TernaryPropertyContainer<double>*>	ternary_materials;
		
};

} // namespace negf

#endif /*MATERIALDATABASE_H_NEGF*/
