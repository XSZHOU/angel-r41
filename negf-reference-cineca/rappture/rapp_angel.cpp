/* compile using:
g++ -I ../../rappture/20090818/include \
	-o rapp_angel.bin rapp_angel.cpp \
	-DNEGF_DIR="\/home/steiger/etc\" \
	-Bdynamic -Wl,-rpath,../../rappture/20090818/lib -L ../../rappture/20090818/lib -lrappture
*/

/* ========================================================

   ANGEL - A NEGF Simulator aimed at LEDs
   
   This file provides a wrapper script to be called from 
   a rappture tool.xml file.
   It does the following:
   1. Transform the contents of driver*.xml, created by the user upon 
      clicking the "Simulate" button, into an ANGEL input deck
   2. Run ANGEL using that input deck and store everything in the 
      same directory
   3. Call save_pngs.py, which lets matplotlib generate pictures
      of contour plots (since rappture doesn't have such a widget.
      Do this for all computed voltages
   4. Hook up the generated images in <sequence> sections within
      the <output> section.
   5. Save run*.xml, which will be visualized by Rappture.
   6. Clean up.
   
   Important note: the paths to angel.bin and save_pngs.py need to
   be hard coded because the file is executed from a different 
   location.
   
   ============================================================
*/

#include "rappture.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
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


/*
// for parallel jobs - NYI
void create_submit_script() {
	ofstream fout("submit_script.sh"); // PWD....
	fout << "#!/bin/sh\n\n";
	fout << "trap cleanup HUP INT QUIT ABRT TERM\n\n";
	fout << "cleanup()\n";
	fout <<  "{\n";
	fout << "   kill -TERM `jobs -p`\n";
	fout << "   exit 1\n";
	fout << "}\n\n";   
	
	fout << "cd [pwd]\n";
	fout << "submit -v cluster -n $nodes -w $walltime\ \n";
	fout << "       COMMAND ARGUMENTS &\n";
	fout << "sleep 5\n";
	fout << "wait\n";
	fout.close();

	system("chmod 755 submit_script.sh");
}
*/


int main(int argc, char ** argv) 
{
	rpUtilsProgress( 0,"Starting Rappturized ANGEL...\n");
	
	for (int ii=0; ii<argc; ii++) {
		msg << "argument " << ii << ": " << argv[ii] << "\n";
	}
	ASSERT(argc==2);
	string driver(argv[1]); // last 4 chars are ".xml"
	//driver = "nanohub";
	driver = driver.substr(0, driver.length()-4);
	msg << "Driver: " << driver << "\n";
	

    char buf[1000];

	// -------------------------------
	// get rappture data as an object
	// -------------------------------
    RpLibrary  * lib = rpLibrary(argv[1]);
	ASSERT(lib!=NULL);
	
	// ------------------------------------------------
	// get input data and write simulator input deck
	// ------------------------------------------------
	
	rpUtilsProgress( 1,"Getting user options...\n");
	
    double T=300; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).number(temperature).current", &T));
	
    double Emin=1; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(numerics).number(Emin).current", &Emin));
    double Emax=2; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(numerics).number(Emax).current", &Emax));
    double NE=200; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(numerics).integer(NE).current", &NE));
	
    double V0=0;    ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(voltages).number(V0).current", &V0));
    double V1=0.1;  ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(voltages).number(V1).current", &V1));
    double dV=0.05; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(params).group(voltages).number(dV).current", &dV));
	
    double l0_mat=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).choice(mat0).current"   , &l0_mat));
    double l0_len=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(length0).current", &l0_len));
    double l0_dx=0;  ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(dx0).current"    , &l0_dx));
	double l0_dop=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(doping0).current", &l0_dop));
	
    double l1_mat=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).choice(mat1).current"   , &l1_mat));
    double l1_len=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(length1).current", &l1_len));
    double l1_dx=0;  ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(dx1).current"    , &l1_dx));
	double l1_dop=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(doping1).current", &l1_dop));
	
    double l2_mat=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).choice(mat2).current"   , &l2_mat));
    double l2_len=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(length2).current", &l2_len));
    double l2_dx=0;  ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(dx2).current"    , &l2_dx));
	double l2_dop=0; ASSERT(0==rpGetDouble(lib,"input.group(tabs).group(struct).structure.current.parameters.number(doping2).current", &l2_dop));
	
	rpUtilsProgress( 2, "Creating input deck...\n");
	ofstream cmd;
	//cmd.open("nanohub.cmd");
	string s = driver+".cmd";  cmd.open(s.c_str());
	
	cmd << "options { \n";
	cmd << "\n";
	cmd << "    temperature                   = " << T << "\n";
	cmd << "\n";
	cmd << "	kp_method                     = 0  \n";
	cmd << "	maximal_k_value               = 1.3 # in 1/nm \n";
	cmd << "	num_k_points                  = 1  \n";
	cmd << "\n";
	cmd << "	num_energy_points             = " << NE << " \n";
	cmd << "\n";
	cmd << "	min_energy                    = " << Emin << " \n";
	cmd << "	max_energy                    = " << Emax << " \n";
	cmd << "\n";
	cmd << "	StrainPolarization 			  = 0  \n";
	cmd << "  \n";
	cmd << "	PotentialUnderrelaxation      = 0.1  \n";
	cmd << "	PulayMixing                   = 0  \n";
	cmd << "	Resonances                    = 0  \n";
	cmd << "	PMLResonances                 = 0  \n";
	cmd << "	IncludeImaginaryContactStates = 0  \n";
	cmd << "	IncoherentContacts            = 1  \n";
	cmd << "	IncoherentContactBroadening   = 1e-2  \n";
	cmd << "\n";
	cmd << "	# ---------------  \n";
	cmd << "	# Self Energies    \n";
	cmd << "	# ---------------  \n";
	cmd << "	ScatteringDecreaseFactor    = 0.05	 # S(0) = (this)*S(end)  \n";
	cmd << "	ScatteringRampFactor        = 2.0	 # S(i) * (this)*S(i-1)  \n";
	cmd << "\n";
	cmd << "	Buettiker                   = 0  \n";
	cmd << "	BuettikerParameter          = 0.01   # [eV]  \n";
	cmd << "	GolizadehMomentumRelaxation = 0  \n";
	cmd << "	GolizadehMomentumParameter  = 0.0001 # [eV^2]  \n";
	cmd << "	GolizadehDephasing          = 0  \n";
	cmd << "	GolizadehDephasingParameter = 0.0001 # [eV^2] sqrt(0.0001)=10meV  \n";
	cmd << "\n";
	cmd << "	OpticalPhonons              = 0 \n";
	cmd << "	LuisierSRpop                = 0 \n";
	cmd << "\n";
	cmd << "	AcousticPhonons             = 0 \n";
	cmd << "\n";
	cmd << "	IonizedImpurities           = 0 \n";
	cmd << "\n";
	cmd << "	# numerical options \n";
	cmd << "	inner_errcrit               = 1e-6 \n";
	cmd << "} \n";
	cmd << "\n";
	cmd << "\n";
	cmd << "experiment_0 { \n";
	cmd << "	Source_voltage = 0.0 \n";
	cmd << "	Drain_min      = -" << V0 << " \n";
	cmd << "	Drain_max      = -" << V1+0.001 << " \n";
	cmd << "	Drain_step     = -" << dV << " \n";
	cmd << "	PotentialUnderrelaxation = 0.1 \n";
	cmd << "	LateUnderrelaxation      = 0.0 \n";
	cmd << "} \n";
	cmd << " \n";
	cmd << " \n";
	cmd << "regions { \n";
	cmd << "	region0_length = " << l0_len << " \n";
	cmd << "	region0_dx     = " << l0_dx  << " \n";
	cmd << "	region0_mat    = " << l0_mat << " \n";
	cmd << "	region0_molefr = 0.3 \n";
	cmd << "	region0_doping = " << l0_dop << " \n";
	cmd << "	\n";
	cmd << "	region1_length = " << l1_len << " \n";
	cmd << "	region1_dx     = " << l1_dx  << " \n";
	cmd << "	region1_mat    = " << l1_mat << " \n";
	cmd << "	region1_molefr = 0.3 \n";
	cmd << "	region1_doping = " << l1_dop << " \n";
	cmd << "	 \n";
	cmd << "	region2_length = " << l2_len << " \n";
	cmd << "	region2_dx     = " << l2_dx  << " \n";
	cmd << "	region2_mat    = " << l2_mat << " \n";
	cmd << "	region2_molefr = 0.3 \n";
	cmd << "	region2_doping = " << l2_dop << " \n";
	cmd << "} \n";
		
	cmd.close();
	
	// ------------------------
	// run simulator
	// ------------------------
	rpUtilsProgress( 5, "Running ANGEL...\n"); 

	msg << "Running simulation...";
	//system("../angel.bin nanohub angel.log ."); // 2nd argument = logfile, 3rd argument = output directory
	s = NEGF_DIR; s.append("/bin/angel.bin  "+driver+"  "+driver+".log  ./"); 
	msg << "ANGEL call: " << s; 
	system(s.c_str());
	msg << " done.\n";
		
	// ------------------------------------------------
	// get simulation output and put back into rappture
	// ------------------------------------------------
	rpUtilsProgress(80, "Retrieving simulation output...\n");
	if (V1!=V0) {
	   ASSERT(fabs(dV)>1e-9 && sign(dV)==sign(V1-V0));
	}
	
	int num_steps = 0;
	vector<double> voltages;
	double volt;
	for(volt=V0; fabs(volt-V0)<=fabs(V1-V0); volt+=dV)
	{
		num_steps++;
		msg << "V0=" << V0 << ", dV=" << dV << ", V1=" << V1 << ", volt=" << volt << "\n";
		voltages.push_back(volt);
		
	    msg << "processing voltage " << volt << "\n";

        msg << "image processing using MATLAB...\n";
        
        //system("matlab -nodesktop -nosplash -r save_pngs"); // saves .png's in calling directory
        //system("/apps/rhel5/MATLAB_R2010a/bin/matlab -nodisplay -nosplash -r save_pngs");
        //system("matlab -nosplash -r save_pngs");
        //system("python save_pngs.py"); // matplotlib
		sprintf(buf, "%s_Drain-%.3fV", driver.c_str(), volt);
		s = "python "; s.append(NEGF_DIR); s.append("/rappture/save_pngs.py ./ "); s.append(buf);
        system(s.c_str()); // matplotlib

        msg << "transfering images to RAPPTURE...\n";
		
		char buf2[1000];
		
        //sprintf(buf,"%d", num_steps);
        sprintf(buf,"%.3f", volt);
		sprintf(buf2, "output.sequence(nE).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));
		sprintf(buf2, "output.sequence(LDOS).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));
		sprintf(buf2, "output.sequence(JE).element(%d).index", num_steps);
		ASSERT(0==rpPutString(lib, buf2, buf, 0));

		sprintf(buf,"%s_Drain-%.3fV_npE.png", driver.c_str(), volt); 
		sprintf(buf2, "output.sequence(nE).element(%d).image(img_nE).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));  // compress, append
		sprintf(buf2, "output.sequence(nE).element(%d).image(img_nE).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "n(x,E)", 0));

        sprintf(buf,"%s_Drain-%.3fV_LDOS.png", driver.c_str(), volt); 
		sprintf(buf2, "output.sequence(LDOS).element(%d).image(img_LDOS).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));
		sprintf(buf2, "output.sequence(LDOS).element(%d).image(img_LDOS).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "LDOS(x,E)", 0));

		sprintf(buf,"%s_Drain-%.3fV_JE.png", driver.c_str(), volt); 
		sprintf(buf2, "output.sequence(JE).element(%d).image(img_JE).current", num_steps);
        ASSERT(0==rpPutFile  (lib, buf2, buf, 1,       0));
		sprintf(buf2, "output.sequence(JE).element(%d).image(img_JE).about.label", num_steps);
        ASSERT(0==rpPutString(lib, buf2, "J(x,E)", 0));
    }
	volt -= dV; // now we have last computed voltage stored in volt
	
	// ------------
	// I-V plot
	// ------------
	sprintf(buf, "%s_Drain-%.3fV.ccurrent", driver.c_str(), volt);
	ofstream t("test.txt"); t << buf; t.close();
	
	ifstream iv_file; iv_file.open(buf); ASSERT(iv_file);
	// first 2 lines are comments
	iv_file.getline(buf, 1000); msg << buf << "\n";
	iv_file.getline(buf, 1000); msg << buf << "\n";
	
	// get num_steps lines
	// use 4'th column of each line
	double dump;
	vector<double> current(num_steps, 0.0);
	for (int ii=0; ii<num_steps; ii++) {
		iv_file >> dump; iv_file >> dump; iv_file >> dump; iv_file >> current[ii]; iv_file >> dump;
	}
	iv_file.close();
	
	msg << "IV-data: ";
	for (uint ii=0; ii<current.size(); ii++) msg << current[ii] << "  ";
	msg << "\n";
	
	// put into rappture
	ostringstream iv_stream;
	for (int ii=0; ii<num_steps; ii++) {
        iv_stream << voltages[ii] << "  " << current[ii] << "    ";
	}
    ASSERT(0==rpPutString(lib, "output.curve(IV).component.xy", iv_stream.str().c_str(), 0));
	
	// --------------
	// log file
	// --------------
	
	/*
	string log;	string line;
	ifstream log_stream("angel.log");
	while(std::getline(log_stream,line)) log += line;
	ASSERT(0==rpPutString(lib, "output.log", log.c_str(), 0));
	*/
		
	// ------------------------
	// save to .xml-file
	// ------------------------
    msg << "lib = " << lib << "\n";	
	const char * buf2 = NULL;
	ASSERT(0==rpXml(lib,&buf2));
    msg.outfile << "Created XML file:\n-----------------\n" << buf2 << "\n";
	
	msg << "saving to run???.xml...\n";
	ASSERT(0==rpResult(lib));

	// ------------------------
	// clean up
	// ------------------------
	s = "rm ./"+driver+"_step*";
	system(s.c_str());
	// maybe remove all driver* files?
    ASSERT(0==rpFreeLibrary(&lib));
	return 0;
}

