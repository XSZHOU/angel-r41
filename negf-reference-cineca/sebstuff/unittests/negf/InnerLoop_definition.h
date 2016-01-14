#ifndef FLENS_DEFINITION_H_
#define FLENS_DEFINITION_H_

// ------------------------------------------
// standard includes
// ------------------------------------------
#include "all.h"
#include "cxxtest/TestSuite.h"

// ------------------------------------------
// class includes
// ------------------------------------------
// flens.h is included in all.h
typedef flens::GeMatrix<flens::FullStorage<cplx, flens::ColMajor> > 	CGEMatrix; 
typedef flens::DenseVector<flens::Array<cplx> > 						CDEVector; 
typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 	DGEMatrix; 
typedef flens::DenseVector<flens::Array<double> > 						DDEVector; 

#include "Hamiltonian.h"
#include "Overlap.h"
#include "NEGFObject.h"
#include "GreenFunctions.h"
#include "SelfEnergies.h"
#include "SEContacts.h"
#include "PostProcessing.h"
#include "InnerLoop.h"

using namespace negf;
using namespace flens;

namespace negf {

class InnerLoopTest : public CxxTest::TestSuite 
{
public:
	
	std::fstream fout;	
	void setUp() {
		N = 20;
		M = 30;
		eps = 1e-13;	
	}
	void tearDown() {
	}
		
	void test_InnerLoop()
	{
		fout.open("./negf/InnerLoop_output.log", std::ios::out); 
		// current directory is Makefile-directory
		logmsg->add_listener(fout);
		
		// test calculation of GF
		void calculate_retarded_gf();
		void calculate_lesser_gf();
		void calculate_greater_gf();
		fout.close();	
	}
	
};

} // end namespace negf

#endif /*FLENS_DEFINITION_H_*/
