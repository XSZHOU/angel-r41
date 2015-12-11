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
#ifndef TERNARYPROPERTYCONTAINER_H_NEGF
#define TERNARYPROPERTYCONTAINER_H_NEGF

#include "all.h"

#include "PropertyContainer.h"

using namespace std;

namespace negf {

	/** property container for ternary materials (AlGaAs etc.) */
	template<class T>
	class TernaryPropertyContainer: public PropertyContainer<T>
	{
	public:
		TernaryPropertyContainer(const char* filename, const double molefraction_);
		virtual ~TernaryPropertyContainer();

		double get_molefraction() const { return molefraction; }
		const PropertyContainer<T> * get_first_pure_material()  const { return first_pure_material; }
		const PropertyContainer<T> * get_second_pure_material() const { return second_pure_material; }

		// ------------------------------------------------------
		// access functions overwritten from PropertyContainer
		// ------------------------------------------------------
		virtual bool valid() 								 const throw(Exception*);
		virtual bool valid_value(const char* key,   const T & value) const throw(Exception*); 
		virtual bool valid_value(const string &key, const T & value) const throw(Exception*); 
		virtual bool valid_key	(const char* key)			 const throw(Exception*);
		virtual bool valid_key	(const string& key)			 const throw(Exception*);
		virtual bool is_set		(const char* key)			 const throw(Exception*);
		virtual bool is_set		(const string& key)			 const throw(Exception*);	
		virtual T  	 get		(const char* key)			 const throw(Exception*);
		virtual T 	 get		(const string& key)			 const throw(Exception*);
		

	protected:
		double molefraction;
		
		PropertyContainer<T> * first_pure_material;		// x=0
		PropertyContainer<T> * second_pure_material;	// x=1
		PropertyContainer<T> * molefraction_dependence;	// additional stuff apart from linear interpolation between x=0,1
	};


template<class T>
TernaryPropertyContainer<T>::TernaryPropertyContainer(const char* filename_char, const double molefraction_):
	molefraction(molefraction_)
{STACK_TRACE(
	NEGF_ASSERT(molefraction>=0.0 && molefraction <= 1.0, "molefraction must be between 0 and 1.");
	
	// divide input filename into bare file and directory
	string filename(filename_char);
	string barefile(filename);
	int loc = -1;
	do
	{
		barefile = barefile.substr(loc+1, barefile.size());
		loc = barefile.find("/", 0);
	} while (loc != (int) string::npos);
	NEGF_ASSERT(barefile.size()>0, "Invalid filename.");
	string filedir = filename.substr(0,filename.size()-barefile.size());
	
	// strip filename to its basics
	string first_pure_name("");
	string second_pure_name("");
	string alloy_name("");
	NEGF_ASSERT(   Constants.ternary_names.size()==Constants.first_pure_names.size()
				&& Constants.ternary_names.size()==Constants.second_pure_names.size(), 
				"check ternary & pure names lists!");
	for (uint ii = 0; ii < Constants.ternary_names.size(); ii++)
	{
		 if (barefile.find(Constants.ternary_names[ii]) != string::npos) { 
		 	first_pure_name  = filedir + Constants.first_pure_names[ii]  + ".mat"; 
		 	second_pure_name = filedir + Constants.second_pure_names[ii] + ".mat"; 
		 	alloy_name       = filedir + Constants.ternary_names[ii]     + ".mat"; 
		 }
	}
	NEGF_FASSERT( first_pure_name!="" && second_pure_name!="" && alloy_name!="",
			"%s is not recognized as a ternary material by TernaryPropertyContainer.", barefile.c_str() );
	NEGF_FASSERT(alloy_name.compare(filename)==0, " %s != %s", alloy_name.c_str(), filename_char);
	
	this->first_pure_material     = new PropertyContainer<T>();
	first_pure_material->read_data_from_file(first_pure_name.c_str());
	this->second_pure_material    = new PropertyContainer<T>();
	second_pure_material->read_data_from_file(second_pure_name.c_str());
	this->molefraction_dependence = new PropertyContainer<T>();
	molefraction_dependence->read_data_from_file(filename_char);
);}

template<class T>
TernaryPropertyContainer<T>::~TernaryPropertyContainer()
{STACK_TRACE(
	delete this->first_pure_material;
	delete this->second_pure_material;
	delete this->molefraction_dependence;
);}
	

/** returns true when all mandatory fields are set in both constituting binary materials */
template<class T>
bool TernaryPropertyContainer<T>::valid() const throw(Exception*)
{STACK_TRACE(
	return (first_pure_material->valid() && second_pure_material->valid());
);}

template<class T>
bool TernaryPropertyContainer<T>::valid_value(const char* key, const T & value) const throw(Exception*)
{STACK_TRACE(
	string tmp(key);
	return this->valid_value(tmp, value);
);}

/** returns true when the value is inside the error/warning bounds of both two constituting binary materials */
template<class T>
bool TernaryPropertyContainer<T>::valid_value(const string &key, const T & value) const throw(Exception*)
{STACK_TRACE(
	if (!(first_pure_material->valid_value(key, value) && second_pure_material->valid_value(key, value))) {
		return false;
	} else { 
		return true;
	}
);}

template<class T>
bool TernaryPropertyContainer<T>::valid_key(const char* key) const throw(Exception *)
{
	string tmp(key);
	return this->valid_key(tmp);	
}

/** returns true when the key is set for both two constituting binary materials since a value
 *  can be calculated if it is set in that case */
template<class T>
bool TernaryPropertyContainer<T>::valid_key(const string& key) const throw(Exception *)
{
	/* hack for bandgap bowing */
	if (key.compare("bandgap_bowing")==0) {
		return molefraction_dependence->valid_key("bandgap_bowing");
	}
	if (key.compare("bandgap_bowing_x")==0) {
		return molefraction_dependence->valid_key("bandgap_bowing_x");
	}
	
	return first_pure_material->valid_key(key) && second_pure_material->valid_key(key);
}

template<class T>
bool TernaryPropertyContainer<T>::is_set(const char* key) const throw(Exception*) 
{
	string tmp(key);
	return this->is_set(tmp);
}

/** returns true when the key is set for both two constituting binary materials since a value
 *  can be calculated in that case */
template<class T>
bool TernaryPropertyContainer<T>::is_set(const string& key) const throw(Exception*) 
{
	/* hack for bandgap bowing */
	if (key.compare("bandgap_bowing")==0) {
		return molefraction_dependence->is_set("bandgap_bowing");
	}
	if (key.compare("bandgap_bowing_x")==0) {
		return molefraction_dependence->is_set("bandgap_bowing_x");
	}
	
	return first_pure_material->is_set(key) && second_pure_material->is_set(key);
}	

template<class T>
T TernaryPropertyContainer<T>::get(const char* key) const throw(Exception*) 
{
	string tmp(key);
	return this->get(tmp);
}

/** returns linearly interpolated value :
 *  A(x) = (1-x)*A(0) + x*A(1) 
 *  ... EXCEPT when A_x2 is set (bowing parameter), in which case
 *  A(x) = (1-x)*A(0) + x*A(1) - x*(1-x)*A_x2
 * 
 *  NOTE that the bowing parameter is the same when x is replaced by 1-x
 *  NOTE the minus sign!!!
 */
template<class T>
T TernaryPropertyContainer<T>::get(const string& key) const throw(Exception*) 
{
	double x = this->molefraction;
	
	/* hack for bandgap bowing */
	if (key.compare("bandgap_bowing")==0) {
		NEGF_FASSERT(molefraction_dependence->is_set("bandgap_bowing"), "key \"bandgap_bowing\" was not found in material %s.", 
						molefraction_dependence->get_name().c_str() );
		double result = molefraction_dependence->get("bandgap_bowing");
		
		// bowing parameter itself can be molefraction-dependent
		if (molefraction_dependence->is_set("bandgap_bowing_x")) {
			//cout << "result = " << result << " += " << molefraction * molefraction_dependence->get("bandgap_bowing_x") << endl;
			result += x * molefraction_dependence->get("bandgap_bowing_x");
		}
		return result;
	}
	
	
	NEGF_FASSERT(first_pure_material->is_set(key), "key \"%s\" was not found in material %s.",key.c_str(), first_pure_material->get_name().c_str());
	NEGF_FASSERT(second_pure_material->is_set(key), "key \"%s\" was not found in material %s.",key.c_str(), second_pure_material->get_name().c_str());
		
	double value_0 = first_pure_material->get(key);
	double value_1 = second_pure_material->get(key);
	string key_x2 = key + "_x2";
	
	/* hacks for Arora model, see Sootodeh et al., J. Appl.Phys. 87, 2890 (2000) */
	if (key.compare("mobility_arora_AN_e")==0 || key.compare("mobility_arora_AN_h")==0) {
		NEGF_ASSERT(!molefraction_dependence->is_set(key_x2), "bowing of this parameter is not used.");
		
		// power interpolation
		return negf_math::pow(value_0, 1.0-x) * negf_math::pow(value_1, x);
	}
	if (key.compare("mobility_arora_alphad_e")==0 || key.compare("mobility_arora_alphad_h")==0) {
		NEGF_ASSERT(!molefraction_dependence->is_set(key_x2), "bowing of this parameter is not used.");
		
		double frac = 1.0 / (1.0 + x * (1.0-x));	// =1 for x=0,1
		return frac * ((1.0-x)*value_0 + x*value_1);
	}
	if (this->get_name().substr(0,6)=="AlGaAs" && 
		(key.compare("mobility_arora_Amin_e")==0 || key.compare("mobility_arora_Ad_e")==0)) {
		NEGF_ASSERT(!molefraction_dependence->is_set(key_x2), "bowing of this parameter is not used.");
		
		double me_0 		= first_pure_material ->get("electron_effective_mass");
		double me_x 		= this				  ->get("electron_effective_mass");
		double eps_stat_0 	= first_pure_material ->get("static_dielectric_constant");
		double eps_stat_x 	= this				  ->get("static_dielectric_constant");
		double eps_opt_0 	= first_pure_material ->get("optic_dielectric_constant");
		double eps_opt_x 	= this				  ->get("optic_dielectric_constant");
		
		double me_factor = negf_math::pow(me_0/me_x, 1.5);
		double eps_factor = (1.0/eps_opt_0 - 1.0/eps_stat_0) / (1.0/eps_opt_x - 1.0/eps_stat_x);
		
		//cout << "value_0 = " << value_0 << ", value_x=" << value_0 * me_factor * eps_factor << endl;
		return value_0 * me_factor * eps_factor;
	}
	
	double alloy_value = (1.0-molefraction)*value_0 + molefraction*value_1;
	
	if (molefraction_dependence->is_set(key_x2)) {
		double bowing = molefraction_dependence->get(key_x2);
		alloy_value -= molefraction * (1.0-molefraction) * bowing;
	}
	
	return alloy_value;
}

	
} // end of namespace negf

#endif /*TERNARYPROPERTYCONTAINER_H_NEGF*/
