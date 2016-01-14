/* compile using:
g++ -I ../../rappture/20090818/include \
	-o see_results.bin see_results.cpp \
	-DNEGF_DIR="\/home/steiger/etc\" \
	-Bdynamic -Wl,-rpath,$PWD/../../rappture/20090818/lib -L $PWD/../../rappture/20090818/lib -lrappture
	
paths to see_results.xml and save_pngs.py can be hardcoded this way, and the binary can 
be called from another place (e.g. the place containing the simulation results)
*/


/*
  This is a tool to visualize ANGEL results using the Rappture GUI even
  when the results were not produced within a GUI.
  
  The simulation to visualize is handed over as the first command line argument.
  
  The tool looks for the computed voltages in the output directory.
  It then calls matplotlib to generate images and hooks them up with the Rappture widgets. 
  The resulting .xml-file is saved in visualize.xml using some dummy .xml-driver file as 
  template (see_results.xml).
  
  In the end a system call starts up rappture using the "rerun" command with visualize.xml
  as an argument.
  
  Example: when in the folder with files nano_Source0.000V_nE, nano_Source0.050V_JE etc., call
  
  ../../rappture/see_results.bin nano

  Copyright (c) 2010, Sebastian Steiger, steiger@purdue.edu
*/

#include "rappture.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <sys/stat.h>

int sign(double x) { return x>=0; }
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::vector;
using std::string;
//using namespace std;

// Messaging class to redirect output to cout AND a file
class Messager {
public:
	Messager()  { outfile.open("rapp_angel.log"); }
	~Messager() { outfile.close(); }	
	template<class T> Messager & operator<<(T argument) { std::cout << argument; outfile << argument; return *this; }	
	static Messager & endl(Messager & stream) { stream.line_break(); return stream; }		
	ofstream outfile; // file keeps track of what is happening here	
private:
	void line_break() { std::cout << std::endl; outfile << std::endl; }
};
Messager msg;


#define ASSERT(x) if (!(x)) { msg << "assertion failed: " << #x << Messager::endl; exit(1); }



int main(int argc, char ** argv) 
{	
	// hand over output-file collection as argument
	for (int ii=0; ii<argc; ii++) {
		msg << "argument " << ii << ": " << argv[ii] << "\n";
	}
	ASSERT(argc==2);	
	string driver(argv[1]);
	msg << "Visualize: " << driver << "\n";
	

    char buf[1000];
	string s;

	// -------------------------------
	// get rappture data as an object
	// -------------------------------
	s = NEGF_DIR; s.append("/rappture/see_results.xml");
    RpLibrary  * lib = rpLibrary( s.c_str() );
	ASSERT(lib!=NULL);
	
	// ------------------------------------------------
	// check which results are present
	// ------------------------------------------------
	
	int num_steps = 0;
	vector<double> voltages;
	vector<string> voltage_names;
	double volt;
	struct stat tmpstat;
	
	// check for voltages between -2.000V, -1.999V, ..., 1.999V, 2.000V
	bool converged = false;
	for (volt=-2.000; volt<2.001; volt+=0.001) {
		if (fabs(volt)<1e-6) { volt = +1.0; volt = 0.0; } // correct 0.0 sign
		sprintf(buf, "%s_Drain%.3fV_nE", driver.c_str(), volt);
		
		if (stat(buf, &tmpstat)==0) {
			voltages.push_back(volt); 
			sprintf(buf, "%s_Drain%.3fV", driver.c_str(), volt);
			voltage_names.push_back(buf); 
			converged = true;
		}
	}
	for (volt=-2.000; volt<2.001; volt+=0.001) {
		if (fabs(volt)<1e-6) { volt = +1.0; volt = 0.0; } // correct 0.0 sign
		sprintf(buf, "%s_Source%.3fV_nE", driver.c_str(), volt);
		
		if (stat(buf, &tmpstat)==0) {
			voltages.push_back(volt); 
			sprintf(buf, "%s_Source%.3fV", driver.c_str(), volt);
			voltage_names.push_back(buf); 
			converged = true;
		}
	}
	
	// if there is no converged voltage, look for <name>_NOTCONVERGED
	if (voltages.size()==0) {
		sprintf(buf, "%s_NOTCONVERGED_nE", driver.c_str(), volt);
		if (stat(buf, &tmpstat)==0) {
			voltages.push_back(-888.888); 
			sprintf(buf, "%s_NOTCONVERGED", driver.c_str(), volt);
			voltage_names.push_back(buf); 			
		}
	}
	
	// if we still can't find anything, look for <name>_stepN
	if (voltages.size()==0) {
		for (uint ii=0; ii<30; ii++) {
			sprintf(buf, "%s_step%d_nE", driver.c_str(), ii);
			if (stat(buf, &tmpstat)==0) {
				voltages.push_back(ii); 
				sprintf(buf, "%s_step%d", driver.c_str(), ii);
				voltage_names.push_back(buf); 			
			}
		}
	}
	
	if (voltages.size()==0) {
		msg << "no result files could be found. aborting.\n";
		exit(1);
	}
	
	
	// ------------------------------------------------
	// get simulation output and put back into rappture
	// ------------------------------------------------
	
	for(uint vv=0; vv<voltages.size(); vv++)
	{
		volt = voltages[vv];
		const char * bufV = voltage_names[vv].c_str();
	
		num_steps++;
		
	    msg << "processing voltage " << volt << "\n";

        msg << "image processing using MATLAB...\n";
		
		
		s = "python "; s.append(NEGF_DIR); s.append("/rappture/save_pngs.py ./ "); s.append(bufV);
        system(s.c_str()); // matplotlib

        msg << "transfering images to RAPPTURE...\n";
		
		char buf2[1000];
		
        sprintf(buf,"%.3f", volt);
		sprintf(buf2, "output.sequence(nE).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));
		sprintf(buf2, "output.sequence(LDOS).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));
		sprintf(buf2, "output.sequence(JE).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));

		sprintf(buf,"%s_npE.png", bufV, volt); 
		sprintf(buf2, "output.sequence(nE).element(%d).image(img_nE).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));  // compress, append
		sprintf(buf2, "output.sequence(nE).element(%d).image(img_nE).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "n(x,E)", 0));
 
        sprintf(buf,"%s_LDOS.png", bufV, volt); 
		sprintf(buf2, "output.sequence(LDOS).element(%d).image(img_LDOS).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));
		sprintf(buf2, "output.sequence(LDOS).element(%d).image(img_LDOS).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "LDOS(x,E)", 0));
 
		sprintf(buf,"%s_JE.png", bufV, volt); 
		sprintf(buf2, "output.sequence(JE).element(%d).image(img_JE).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));
		sprintf(buf2, "output.sequence(JE).element(%d).image(img_JE).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "J(x,E)", 0));
    }
	
	// ------------
	// I-V plot
	// ------------
	if (converged) 
	{
		// take file with maximum ABSOLUTE value
		uint maxidx = 0;
		double maxvolt = 0.0;
		for (uint ii=0; ii<voltages.size(); ii++) {
			if (fabs(voltages[ii])>maxvolt) {
				maxidx = ii;
				maxvolt = fabs(voltages[ii]);
			}
		}
		
		
		sprintf(buf, "%s.ccurrent", voltage_names[maxidx].c_str()); // last file
	
		ifstream iv_file; iv_file.open(buf); ASSERT(iv_file);
		// first 2 lines are comments
		iv_file.getline(buf, 1000); msg << buf << "\n";
		iv_file.getline(buf, 1000); msg << buf << "\n";
	
		// get num_steps lines
		// use 4'th column of each line
		double dump;
		vector<double> iv_volt(num_steps, 0.0);
		vector<double> current(num_steps, 0.0);
		for (int ii=0; ii<num_steps; ii++) {
			iv_file >> iv_volt[ii]; iv_file >> dump; iv_file >> dump; iv_file >> current[ii]; iv_file >> dump;
		}
		iv_file.close();
	
		msg << "IV-data: ";
		for (uint ii=0; ii<current.size(); ii++) msg << current[ii] << "  ";
		msg << "\n";
	
		// put into rappture
		ostringstream iv_stream;
		for (int ii=0; ii<num_steps; ii++) {
    	    iv_stream << iv_volt[ii] << "  " << current[ii] << "    ";
		}
    	ASSERT(0==rpPutString(lib, "output.curve(IV).component.xy", iv_stream.str().c_str(), 0));
	}
	
	// --------------
	// log file
	// --------------
	
		
	// ------------------------
	// save to .xml-file
	// ------------------------
    msg << "lib = " << lib << "\n";	
	const char * buf2 = NULL;
	ASSERT(0==rpXml(lib,&buf2));
    msg.outfile << "Created XML file:\n-----------------\n" << buf2 << "\n";
	
	//msg << "saving to run???.xml...\n";
	//ASSERT(0==rpResult(lib));

	// ------------------------
	// clean up
	// do NOT remove output files
	// ------------------------
    ASSERT(0==rpFreeLibrary(&lib));
	
	// ------------------------
	// run rappture
	// ------------------------
	msg << "saving visualize.xml...\n";
	ofstream fout("visualize.xml");
	fout << buf2;
	fout.close();
	
	msg << "starting rappture GUI... close GUI to exit.\n";
	s = "rerun visualize.xml";
	system(s.c_str());
	return 0;
}

