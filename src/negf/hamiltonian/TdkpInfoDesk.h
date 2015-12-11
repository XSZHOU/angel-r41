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
#ifndef TDKPINFODESK_H_
#define TDKPINFODESK_H_

#include "all.h"

//#include "Bandedge.h"
#include "PropertyContainer.h"
#include "TernaryPropertyContainer.h"
#include "MaterialDatabase.h"

#ifndef NOTDKP
#include "tdkp/interface/Interface.h"  	// interface Tdkp-->Negf   -I$(INC_PATH)/...
#endif

// instead of the standard stack trace, which throws an Exception*, throw a simple std::string which can be caught by tdkp
#define INTERFACE_STACK_TRACE(...) \
	try { __VA_ARGS__ } catch(Exception *e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw e->get_reason(); }


using namespace std;

namespace negf {
	
	/** returns 1.0 or 2.0, the latter in case of effective mass calculation for electrons
	 *  used in many KP classes */
#ifndef NOTDKP
	double get_spin_degeneracy(tdkp::Interface * kpresult, quantities::PhysicalQuantity e_or_h);
#endif
	double get_spin_degeneracy(const double & kpmodel, quantities::PhysicalQuantity e_or_h);
	
#ifndef NOTDKP
	/** Interface Negf-->Tdkp, implements tdkp::InformationDesk */
	class TdkpInfoDesk : public tdkp::InformationDesk
#else
	class TdkpInfoDesk 
#endif
	{
	public:
	
		TdkpInfoDesk(const MaterialDatabase * db_, double temperature_): db(db_), temperature(temperature_) 
			{INTERFACE_STACK_TRACE( NEGF_ASSERT(db!=NULL, "null pointer encountered."); );}
		
	 	virtual ~TdkpInfoDesk()  {INTERFACE_STACK_TRACE();}
	 
		// ---------------------------------------------------------------
		// interface functions
		// ---------------------------------------------------------------
		
	 	// get a material property --> no ambiguities
	 	virtual double get_property(const string & materialname, const string & propertyname) const;
	 	virtual double get_slc_param(const std::string & propertyname) const;
	 	
	 	// check if a property is available
	 	virtual bool is_set(const string & materialname, const string & propertyname) const;
	 		 	
		// special treatment of conduction band edges
		static double get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/, 
									 const MaterialDatabase * const database) /*const*/;	
	 	
	protected:
		// unit converter to TDKP units
		double convert_to_tdkp_units(units::UnitType unit, const double & value) const;
			
	 	const MaterialDatabase * db;
	 	double temperature;	// for band edges
	};
	
} // end of namespace 

#endif /*TDKPINFODESK_H_*/
