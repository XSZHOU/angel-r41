# module name needs to be the same as in _<name>.so !
%module negf

%{
	#include <iostream>
	#include <sys/time.h>
	#include <string>
	#include <sstream>
	#include <fstream>
	#include <stdio.h>
	#include <stdlib.h>
	#include <vector>
	
	#include "ForwardDeclarations.h"
	
	// now include all the files that are %included for SWIG
	
	#include "all.h"
	#include "Logger.h"
	#include "Timer.h"
	#include "Filenames.h"
	#include "Interrupt.h"
	#include "Mat.h"
	#include "MPI.h"
	
	#include "InputParser.h"
	#include "OutputData.h"
	#include "MaterialDatabase.h"
	#include "PropertyContainer.h"
	
	#include "MaterialDatabase.h"
	
	#include "Vertex.h"
	#include "Edge.h"
	#include "Element.h"
	#include "Region.h"
	#include "Contact.h"
	#include "Geometry.h"
	#include "BoxMethod.h"
	
	#include "Equation.h"
	#include "ExplicitEquation.h"
	#include "ImplicitEquation.h"
	#include "ArbitraryData.h"
	#include "Voltages.h"
	#include "Poisson.h"
	
	#include "NewtonSolver.h"
	
	#include "Energies.h"
	#include "Kspace.h"
	
	#include "Options.h"
	#include "NEGFObject.h"
	#include "GreenFunctions.h"
	
	#include "SelfEnergy.h"
	#include "SEContacts.h"
	#include "SEBuettiker.h"
	#include "SEGolizadeh.h"
	#include "SEOpticalPhonon.h"
	#include "SEAcousticPhonon.h"
	#include "SEPhotonSpontaneous.h"
	#include "SEIonizedImpurities.h"
	#include "SelfEnergies.h"
	
	#include "Current.h"
	#include "Luminescence.h"
	#include "PostProcessing.h"
	#include "TdkpInfoDesk.h"
	#include "Hamiltonian.h"
	#include "ContactFermilevel.h"
	#include "Overlap.h"
	#include "InnerLoop.h"
	#include "PoissonProblem.h"
	#include "OuterLoop.h"
%}

# --------------------------
# SWIG library functions
# --------------------------
%include "std_string.i"
%include "std_vector.i"
namespace std {
	%template(IntVector) vector<int>;
	%template(DoubleVector) vector<double>;
	%template(DoubleDoubleVector) vector< vector<double> >;
	%template(StringVector) vector<string>;
	%template(UnitVector) std::vector<negf::units::UnitType>;
#	%template(EquationVector) vector<negf::Equation *>;
#	%template(OutputDataVector) vector<negf::OutputData *>;
#	%template(VertexVector) vector<negf::Vertex *>;
#	%template(RegionVector) vector<negf::Region *>;
}


# --------------------------
# NEGF functions
# --------------------------

%include "ForwardDeclarations.h"

%include "all.h"

# Tell Python what to do with a thrown object
%include "Exception.h"
%typemap(throws) negf::Exception * %{
	PyErr_SetString(PyExc_RuntimeError, $1->get_reason().c_str());
	SWIG_fail;
%}
%typemap(throws) negf::Exception* %{
	PyErr_SetString(PyExc_RuntimeError, $1->get_reason().c_str());
	SWIG_fail;
%}

	
# Python does not know about ostreams --> replace add_listener, del_listener!
# At the moment python also has problems with the overladed emit-functions with varargs
%ignore negf::Logger::add_listener;
%ignore negf::Logger::del_listener;
%naturalvar negf::Logger::emit;
%include "Logger.h"
%extend negf::Logger {
	std::ofstream * add_device(std::string listener)
	{
		if (listener=="std::cout") {
			self->add_listener(std::cout);
			return 0;
		} else {
			std::ofstream * fout = new std::ofstream(listener.c_str(), ios::app);
			if (!(*fout)) {
				self->emit(negf::LOG_ERROR, "logfile %s could not be created/opened.", listener.c_str());
				return 0;
			}
			self->emit(negf::LOG_INFO, "logfile: %s",listener.c_str());
			self->add_listener(*fout);
			return fout;
		}
	}
	void del_device(std::string listener)
	{
		if (listener=="std::cout") {
			self->del_listener(std::cout);
		} else {
			self->emit(negf::LOG_INFO, "if del_device is called with a string, it must be std::cout");
		}
	}
	void del_device(std::ostream * address)
	{
		self->del_listener(*address);
	}
};

# %apply unsigned int { uint32_t }  # should fix uint problem


%include "Timer.h"
%include "Filenames.h"
%include "Interrupt.h"
%ignore negf::Matc::Matc(const GEMatrix::ConstView);
%ignore negf::Matd::Matd(const DGEMatrix::ConstView);
%include "Mat.h"
%include "MPI.h"

%include "InputParser.h"
%extend negf::InputParser {
	std::vector< const std::vector<double> *> get_better_vector(std::vector< std::vector<double> *> & old_vector)
	{
		std::vector< const std::vector<double> *> new_vector;
		for (uint ii=0; ii < old_vector.size(); ii++) {
			new_vector.push_back(old_vector[ii]);
		}
		return new_vector;
	}
	std::vector< const std::vector<double> *> get_best_vector(std::vector< std::vector<double> > & old_vector)
	{
		std::vector< const std::vector<double> *> new_vector;
		for (uint ii=0; ii < old_vector.size(); ii++) {
			new_vector.push_back(&old_vector[ii]);
		}
		return new_vector;
	}
	const std::vector<double> * get_const_ptr(std::vector<double> & vec)
	{
		return &vec;
	}
};
%include "OutputData.h"
%include "MaterialDatabase.h"
%include "PropertyContainer.h"
%template(PropertyContainerDouble) negf::PropertyContainer<double>;

%include "MaterialDatabase.h"

%include "Vertex.h"
%include "Edge.h"
%include "Element.h"
%include "Region.h"
%include "Contact.h"
%include "Geometry.h"
%include "BoxMethod.h"

# In order for Python to recognize that Poisson is derived from Equation, we need to include all 
# classes in between, too --> ExplicitEquation, ImplicitEquation
%include "Equation.h"
%include "ExplicitEquation.h"
%include "ImplicitEquation.h"
%include "ArbitraryData.h"
%include "Poisson.h"
%include "Voltages.h"

%include "NewtonSolver.h"


%include "Energies.h"
%include "Kspace.h"

# must include NEGFObject becuase python needs to know about the base classes!
%include "Options.h"
%include "TdkpInfoDesk.h"
%include "Hamiltonian.h"
%include "Overlap.h"
%include "NEGFObject.h"
%template(NEGFObjectMatc) negf::NEGFObject<negf::Matc>;
%include "GreenFunctions.h"

%include "SelfEnergy.h"
%include "SEContacts.h"
%include "SEBuettiker.h"
%include "SEGolizadeh.h"
%include "SEOpticalPhonon.h"
%include "SEAcousticPhonon.h"
%include "SEPhotonSpontaneous.h"
%include "SEIonizedImpurities.h"
%include "SelfEnergies.h"

%include "Current.h"
%include "Luminescence.h"
%include "PostProcessing.h"
%include "ContactFermilevel.h"
%include "InnerLoop.h"
%include "PoissonProblem.h"
%include "OuterLoop.h"

