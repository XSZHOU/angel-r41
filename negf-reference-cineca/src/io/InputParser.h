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
#ifndef _INPUTPARSER_H_NEGF
#define _INPUTPARSER_H_NEGF

#include "all.h"

#include "Geometry.h"
// includes of external headers are in InputParser.cpp

using namespace std;

namespace negf {
		
	typedef const vector<double> dblvec;
	typedef const vector<const vector<double> *> valsarray;
	
	/** Class that handles input/output. */
	class InputParser{
	public:	
		InputParser() {};
		~InputParser() {};


		Geometry * read_grd(const map< string, PropertyContainer<double> * > * cmdfile);

#ifndef NODFISE
		// ------------------------------------
		// DF-ISE
		// ------------------------------------
		Geometry * 	read_dfise_grd(const char* filename) const throw (Exception *);  // attaches .grd automatically
		void 	   	write_dfise_grd(const char * filename, const Geometry * const grid) const throw (Exception *);
		
		void		read_dfise_dat(const char* filename, const char* gridname,
							const Geometry * const geom, vector<double> & values, 
							const string & fieldname, const string & loc) const throw (Exception *);
		void 		write_dfise_dat( const char* filename, const Geometry * geom, 
							valsarray values_array, 
							const uint num_datasets, 
							const vector<string> & datanames,
							const vector<string> & locations,
							const vector<units::UnitType> & unitypes) const throw (Exception *);
		
		void		read_dfise_plt(const char* filename, vector<string> & fieldnames, 
								vector< vector<double> > & values) throw (Exception *);
		
		void		write_dfise_plt( const char* filename, 
								 const uint num_values, 
								 const uint num_datasets,
								 const vector<string> & datanames,
								 const vector<units::UnitType> & unittypes,
								 const vector<double> & values) const throw (Exception *);
		
		/** check the existence of some field in some .dat-file */
		bool    	check_field_existence(const char* filename, const string & fieldname) const throw (Exception *);
				
		void 		get_fields(const char* filename, vector<string> & fieldnames, 
								vector<string> & locations) const throw (Exception *);

        /** molefraction preparer */
        void    prepare_molefractions(const Geometry * const geom, const string & filename) const;
#endif

		// --------------------------------
		// helpers, other
		// --------------------------------
		
		void		write_xy(const char* filename, 
							 const uint num_values, 
							 const uint num_datasets,
							 const vector<string> & datanames,
							 const vector<units::UnitType> & unittypes,
							 const vector<double> & values) const throw (Exception *);
		
		// these wrappers are necessary because I can't get python to work directly with the all.h routines
		void 		write_matrix(const char* filename, const Matc & matrix) const throw (Exception *) { negf::write_matrix(filename, matrix); }
		void 		write_matrix(const char* filename, const Matd & matrix) const throw (Exception *) { negf::write_matrix(filename, matrix); }
		void 		read_matrix(const char* filename, Matc & matrix) const throw (Exception *) { negf::read_matrix(filename,matrix); }
		void 		write_xE_matrix(const char* filename, const Matd & matrix, const Geometry * xspace, const vector<double> & energies) const throw (Exception *);
		void 		write_current_matrix(const char* filename, const Matd & matrix, const Geometry * xspace, const vector<double> & energies) const throw (Exception *);
		void 		write_phi_n_p_Ec_Ev_J(const char* filename, const Geometry * xspace, 
						const vector<double> & pot,  const vector<double> & edens, const vector<double> & hdens, 
						const vector<double> & curr, const vector<double> & Ec,    const vector<double> & Ev);
		
		/** function to read command file */
		map< string, PropertyContainer<double> * > read_cmd_file(string filename) const throw (Exception *);
		
		/** functions to read position parameters of quantized regions from file */
		void                     strip_whitespaces (string& str) const;
		

	protected:
		/** transform the .cmd-file material number to the material name
		    this is the only place in the code defining the mapping.
			char * name needs to be preallocated. */
		static void number_to_materialname(const int mat_idx, char * name);
	
	
#ifndef NODFISE
		/** helper functions for conversion DF-ISE/VTK element types --> NEGF element types */
		void 	dfise_to_new(const int & dfise_type, element_type::ElementType & new_type, uint & el_dim) const;
		
		/** helper functions for save_grd_file */
		uint    get_dfise_ordered_index(const Edge * edge, const Vertex * first_vertex) const;
#endif
		Edge *  get_edge_containing_specific_vertices(const Element * elem, uint v0, uint v1) const;
		Edge *  get_edge_containing_specific_vertices(const Face    * face, uint v0, uint v1) const;
		
		/** helper functions for molefraction preparer */
		bool 	check_if_ternary_materials_exist(const Geometry * const geom) const;
		vector<double> get_regionwise_constant_xmole(const Geometry * const geom, const vector<double> & xmole_on_vertices) const;
		void 	assign_xmole_to_regions(const Geometry * const geom, const vector<double> & xmole) const;
		
	
	// ----------------------------------------------------
	// standard routines for reading/writing matrices
	// ----------------------------------------------------
	public:
		template<class T>
		T* read_matrix(const char* filename) const throw(Exception *)
		{STACK_TRACE(
			ifstream fin(filename);
			if(fin) {
				logmsg->emit(LOG_INFO_L2, "reading matrix from %s",filename);
				int size_x, size_y;
				fin >> size_x >> size_y;
				T* matrix = new T[size_x * size_y];
				for(int ii = 0; ii < size_x; ii++) {
					for(int jj = 0; jj < size_y; jj++) {
						fin >> matrix[ii * size_y + jj]; 	
					}	
				} 
				int test;
				fin >> test;
				if(test != 88888888) {
					NEGF_EXCEPTION("something fishy happend");	
				}
				fin.close();
				return matrix;
			} else {
				NEGF_FEXCEPTION("file %s could not be opened for read.", filename);
			}
		);}
		
		template<class T>
		void read_matrix(const char* filename, vector< vector<T> > & result) const throw(Exception *)
		{STACK_TRACE(
			ifstream fin(filename);
			result.clear();
			if(fin) {
				logmsg->emit(LOG_INFO_L2, "reading matrix from %s",filename);
				vector<T> emptyline;
				string s;
				//istringstream linestream;
				while (!fin.eof()) {
					result.push_back(emptyline);
					getline(fin,s);
					//linestream.str(s);
					istringstream linestream(s);
					//cout << "string: \""<<linestream.str()<<"\"        ";
					T num;
					vector<T> & ref_result = result[result.size() - 1];
					while (linestream >> num) {
						ref_result.push_back(num);
						//cout << "whoops!!";
					}
					//cout << "line "<<result.size()-1<<": size="<< ref_result.size() << "              ";
				}
				//cout << "result.size()="<<result.size()<<endl;
				//cout << "secondlast row: size()="<<result[result.size()-2].size()<<endl;
				//if (result[result.size()-2].size()>0)
				//	cout << "first entry of secondlast row: " << result[result.size()-2][0] << endl;
				NEGF_ASSERT(result[result.size()-2].size()==1 && result[result.size()-2][0]==88888888,
							"something fishy happend");	
				result.resize(result.size()-2);
				fin.close();
			} else {
				NEGF_FEXCEPTION("file %s could not be opened for read.", filename);
			}
		);}
		
		template<class T>
		void write_matrix(const char* filename, const T* matrix, const int& size_x, const int& size_y) const throw(Exception *)
		{STACK_TRACE(
			ofstream fout(filename);
			if(fout) {
				logmsg->emit(LOG_INFO_L2, "writing matrix to %s",filename);
				fout.precision(12); 
				fout << size_x << "\t" << size_y << "\n";		
				for(int ii = 0; ii < size_x; ii++) {
					for(int jj = 0; jj < size_y; jj++) {
						fout << matrix[ii * size_y + jj] << "\t"; 	
					}	
					fout << "\n";
				} 				
				fout << 88888888;
				fout << "\n";
				fout.close();
			} else {
				NEGF_FEXCEPTION("file %s could not be opened for write.", filename);	
			}
		);}
		
		template<class T>
		void write_matrix(const char* filename, const vector< vector<T> > & matrix) const throw(Exception *)
		{STACK_TRACE(
			ofstream fout(filename);
			if(fout) {
				logmsg->emit(LOG_INFO_L2, "writing matrix to %s",filename);
				fout.precision(12); 
				for(uint ii = 0; ii < matrix.size(); ii++) {
					for(uint jj = 0; jj < matrix[ii].size(); jj++) {
						fout << matrix[ii][jj] << "\t"; 	
					}	
					fout << "\n";
				}
				fout << 88888888 << "\n";
				fout.close();
			} else {
				NEGF_FEXCEPTION("file %s could not be opened for write.", filename);
			}
		);}
	};
	
} // end of namespace

#endif /*_INPUTPARSER_H_NEGF*/
