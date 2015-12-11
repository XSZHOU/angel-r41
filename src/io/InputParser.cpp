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
/*  Note that DF-ISE used to be a proprietary file format used
	in the ISE AG TCAD suite. ISE was bought by Synopsys in 2004 and 
	there is no information about a license, trademark or copyright in
	the current manual. If you have detailed information or think that this class 
	violates any license, trademark or copyright, please contact the 
	author: steiger@purdue.edu. 
*/
#include "InputParser.h"
using namespace negf;


#include <boost/lexical_cast.hpp>	// $(INC_PATH) must have boost included!


/** Read the command file which has simulation parameters and determines he type of experiment,
 *  i.e. which contacts are ramped
 *  The file has ending .cmd and the following format:
 *  parameters {
 *      temperature = 273
 *      ...
 *  }
 *  experiment_0 {
 *  ...
 *  }
 *  experiment_1 {
 *  ...
 *  }
 *  ...
 */
map< string, PropertyContainer<double> * >  InputParser::read_cmd_file(string filename) const throw (Exception *)
{STACK_TRACE(
	ifstream fin;
	filename.append(".cmd");
	fin.open(filename.c_str());
	if(!fin) {
		logmsg->emit(LOG_ERROR,".cmd-file %s not found.",filename.c_str());
		NEGF_FEXCEPTION(".cmd-file %s is necessary for the simulation.",filename.c_str());
		map< string, PropertyContainer<double> * > sim;
		return sim;	
	}
	string str_in;
	string str_stripped;
	string::size_type loc, loc2;
	string item;
	string::iterator it;
	
	logmsg->emit(LOG_INFO_L2,  "reading data from file %s", filename.c_str());
	
	map< string, PropertyContainer<double> * > sim_props;
	// contains PropertyContainers for the parameters as well as for the experiments

	string reffile(fnames->get_materialdirectory()); reffile.append("/cmd.cnf");
	logmsg->emit(LOG_INFO_L2,  "definition-file for identifier properties is %s).",reffile.c_str());
	
	while(!fin.eof()) {
		loc = string::npos;

		// find identifier, ignore comments
		while(loc==string::npos && !fin.eof()) {
			getline(fin,str_in);
			this->strip_whitespaces(str_stripped = str_in);
			if(str_stripped.size() > 0 && str_stripped[0] != '#')
				loc = str_stripped.find("{", 0);	// identifier that marks the beginning of a region spec
		}
		if (fin.eof())
			break;

		// get name and type of information from str_stripped
		string regname = str_stripped.substr(0, loc);
		if (   regname.find(constants::options_name)==string::npos
		    && regname.find("experiment")==string::npos
		    && regname.find("regions")==string::npos) {
			NEGF_FEXCEPTION("Identifier must contain \"%s\" or \"experiment\" or \"regions\".", constants::options_name.c_str());
		}
		NEGF_ASSERT(sim_props.find(regname)==sim_props.end(), "duplicate identifier.");
		logmsg->emit(LOG_INFO_L2,  "adding identifier %s.",regname.c_str());
		sim_props[regname] = new PropertyContainer<double>(reffile.c_str());
		sim_props[regname]->set_name(regname);
		
		// get properties until region-ending identifier appears
		// any properties on the line of the region-ending identifier will be ignored
		while(!fin.eof() && str_stripped.find("}",0)==string::npos) 
		{
			getline(fin,str_in);
			this->strip_whitespaces(str_stripped = str_in);
			loc = str_stripped.find("=",0);			// property identifier
			if (loc==string::npos || str_stripped[0]=='#') {
				continue;
			}
			
			string property = str_stripped.substr(0, loc);
			logmsg->emit(LOG_INFO_L3,  "   ...reading property \"%s\"",property.c_str());
			NEGF_ASSERT(str_stripped.size()>loc,"Theres nothing behind the '='!");
			str_stripped = str_stripped.substr(loc+1);
				
			// parse line (only doubles allowed) until end of line or '#' is reached
			double value;
			try {
				loc2 = str_stripped.find("#");
				item = str_stripped.substr(0, loc2);
				value = boost::lexical_cast<double>(item);
			} catch(boost::bad_lexical_cast &) {
				logmsg->emit(LOG_ERROR,"in file %s: could not parse argument '%s'", filename.c_str(), item.c_str());
				NEGF_EXCEPTION("Parsing error, aborting. Must be double.");
			}
			
			sim_props[regname]->set(property, value);
		}

	} // end while
	fin.close();

	return sim_props;
);}

void InputParser::strip_whitespaces(string& str) const
{STACK_TRACE(
	string tmp;
	for(string::iterator it = str.begin(); it != str.end(); it++) {
		if((*it) != ' ' && (*it) != '\t') {
			tmp.push_back(*it);
		}
	}
	str = tmp;
);}

Edge * InputParser::get_edge_containing_specific_vertices(const Element * elem, uint v0, uint v1) const
{STACK_TRACE(
	NEGF_ASSERT(v0 < elem->get_num_vertices() && v1 < elem->get_num_vertices(), "invalid argument.");
	for (uint ii = 0; ii < elem->get_num_edges(); ii++)
	{
		Edge * edge = elem->get_edge(ii);
		if (   (   edge->get_lower_vertex()==elem->get_vertex(v0) || edge->get_lower_vertex()==elem->get_vertex(v1))
			&& (   edge->get_upper_vertex()==elem->get_vertex(v0) || edge->get_upper_vertex()==elem->get_vertex(v1)) ) 
		{
			return elem->get_edge(ii);
		}
	}
	NEGF_FEXCEPTION("Edge with vertices %d, %d was not found in element %d's edge list.",v0,v1, elem->get_index_global());
	return 0;
);}

Edge * InputParser::get_edge_containing_specific_vertices(const Face * face, uint v0, uint v1) const
{STACK_TRACE(
	NEGF_ASSERT(v0 < face->get_num_vertices() && v1 < face->get_num_vertices(), "invalid argument.");
	for (uint ii = 0; ii < face->get_num_edges(); ii++)
	{
		Edge * edge = face->get_edge(ii);
		if (   (   edge->get_lower_vertex()==face->get_vertex(v0) || edge->get_lower_vertex()==face->get_vertex(v1))
			&& (   edge->get_upper_vertex()==face->get_vertex(v0) || edge->get_upper_vertex()==face->get_vertex(v1)) ) 
		{
			return face->get_edge(ii);
		}
	}
	NEGF_FEXCEPTION("Edge with vertices %d, %d was not found in face %d's edge list.",v0,v1,face->get_index_global());
	return 0;
);}

void InputParser::number_to_materialname(const int mat_type, char * name)
{
	NEGF_ASSERT(name!=NULL, "need to allocate memory.");
    switch (mat_type) {
        case  0: sprintf(name,"%s",   "GaAs"); break;
        case  1: sprintf(name,"%s",   "AlAs"); break;
        case  2: sprintf(name,"%s",   "InAs"); break;
//        case  3: sprintf(name,"%s",    "GaP"); break;
//        case  4: sprintf(name,"%s",    "AlP"); break;
        case  5: sprintf(name,"%s",    "InP"); break;;
        case  6: sprintf(name,"%s",    "GaN"); break;
        case  7: sprintf(name,"%s",    "AlN"); break;
        case  8: sprintf(name,"%s",    "InN"); break;
        case 10: sprintf(name,"%s", "AlGaAs"); break;
        case 11: sprintf(name,"%s", "InGaAs"); break;
        case 12: sprintf(name,"%s", "InAlAs"); break;
//        case 13: sprintf(name,"%s",  "AlGaP"); break;
//        case 14: sprintf(name,"%s",  "InGaP"); break;
//        case 15: sprintf(name,"%s",  "InAlP"); break;
        case 16: sprintf(name,"%s",  "AlGaN"); break;
        case 17: sprintf(name,"%s",  "InGaN"); break;
//        case 18: sprintf(name,"%s",  "InAlN"); break;
        case 20: sprintf(name,"%s","Silicon"); break;
        case 21: sprintf(name,"%s",   "SiO2"); break;
//        case 22: sprintf(name,"%s",     "Ge"); break;
//        case 23: sprintf(name,"%s",   "SiGe"); break;
        default: NEGF_EXCEPTION("material index not recognized.");
    }
}

/** create geometry from .cmd-file:
 *  regions {
 *       name1_mat    = 0     # some integer. 0: GaAs. 1: InAs. 2: ...
 *       name1_length = 4
 *       name1_dx     = 0.2
 *       name1_doping = 1e19
 *       name1_molefr = 0.2   # only relevant for ternary materials
 *   }
 */
Geometry * InputParser::read_grd(const map< string, PropertyContainer<double> * > * cmdfile)
{
    PropertyContainer<double> * regs = 0;
    for (map< string, PropertyContainer<double> * >::const_iterator it = (*cmdfile).begin(); it!=(*cmdfile).end(); it++) {
        if (it->first=="regions") {
            regs = it->second;
            break;
        }
    }
    NEGF_ASSERT(regs!=0, "[InputParser] regions section was not found.");

    // ---------------------------
    // read in region by region
    // ---------------------------

    vector<Vertex *>  vertices;
    int vertex_count = 0;
    vector<Edge *>    edges;
    int edge_count = 0;
    vector<Element *> elements;
    int elem_count = 0;
    vector<Region *> regions;
    int region_count = 0;
    vector<double> molefraction_per_vertex;

    // the two contacts:
    // for the simulation, the first and last region must be called "contact_0" and "contact_1"
    // they must consist of 3 vertices (2 intervals) with the same spacing as the first/last interior element
    // the 3 vertices need to have a location 'e' (as opposed to interior vertices which are 'i' or 'f')
    // location 'e' is the same as internal index -1
    vector<Contact *> contacts;
    contacts.push_back(new Contact("Source"));
    contacts[0]->set_index(0);
    contacts.push_back(new Contact("Drain"));
    contacts[1]->set_index(1);

    double x = 0.0;
    int vl = -1; // stores first vertex within each layer
    char buf[1000];

    // --------------------------------------------------------------------
    // left contact: automatic creation with dx and material from region0
    // - add Region object containing 3 vertices, 2 elements, with spacing = last spacing
    // - add Contact object with only 1 vertex: 3rd-last vertex of structure
    // --------------------------------------------------------------------
    {
        NEGF_ASSERT(regs->is_set("region0_dx") && regs->is_set("region0_mat"), "need a region 0.");
        double dx = constants::convert_from_SI(units::length, 1e-9*regs->get("region0_dx"));
        int mat_type = int(regs->get("region0_mat") + 0.5);
        double molefraction = 0.0; if (regs->is_set("region0_molefr")) molefraction = regs->get("region0_molefr");

        regions.push_back(new Region("contact_0"));
        regions[region_count]->set_index(region_count);
		number_to_materialname(mat_type, buf);
        regions[region_count]->set_material_name(buf);

        for (int vv=0; vv<3; vv++) {
            if (vv==1 || vv==2) x+=dx;

            //cout << "(left contact) vertex_count=" << vertex_count << ": x=" << x << endl;
            vertices.push_back(new Vertex(vertex_count, x));
            vertex_count++;
            vertices[vertex_count-1]->set_index_internal(-1);
            molefraction_per_vertex.push_back(molefraction);

            if (vv==0) continue;

            edges.push_back(new Edge(edge_count, vertices[vertex_count-1], vertices[vertex_count-2]));
            edge_count++;

            elements.push_back(new Element(elem_count, element_type::interval));
            elem_count++;
            elements[elem_count-1]->add_vertex (vertices[vertex_count-1]);
            elements[elem_count-1]->add_vertex (vertices[vertex_count-2]);
            elements[elem_count-1]->add_edge   (edges[edge_count-1]);
            elements[elem_count-1]->set_region (regions[region_count]);
            regions [region_count]->add_element(elements[elem_count-1]);
        }

        contacts[0]->add_vertex(vertices[vertex_count-1]); // single vertex only

        region_count++;
    }

    // --------------------------------------------------------------------
    // device-internal regions
    // --------------------------------------------------------------------

    while (true) {

        // check if there is an option "region<region_counter>_length"
        sprintf(buf, "region%d_length",region_count-1); // -1 because of left contact
        if (!regs->is_set(buf)) {
            break;
        }

        logmsg->emit(LOG_INFO_L1,"reading region %d...", region_count-1);

        // regionN_length is supposed to be in nm
        double length = constants::convert_from_SI(units::length, 1e-9*regs->get(buf));

        // read in dx [nm]
        sprintf(buf, "region%d_dx",region_count-1);
        NEGF_ASSERT(regs->is_set(buf), "could not find dx");
        double dx = constants::convert_from_SI(units::length, 1e-9*regs->get(buf));

        int     num_intervals     = int(length/dx + 0.5);
        double  num_intervals_dbl = length/dx;
        //cout << "num_intervals: int = " << num_intervals << ", double=" << num_intervals_dbl << endl;
        if (fabs(num_intervals - num_intervals_dbl) > constants::convert_from_SI(units::length, 1e-12)) {
            NEGF_EXCEPTION("dx must be fitting into length");
        }

        // read in material type
        sprintf(buf, "region%d_mat",region_count-1);
        NEGF_ASSERT(regs->is_set(buf), "could not find material type");
        double materialtype = regs->get(buf);
        int mat_type = int(materialtype + 0.5);
        //cout << "materialtype=" << materialtype << ", mat_type=" << mat_type << endl;

        // read in material molefraction
        sprintf(buf, "region%d_molefr",region_count-1);
        double molefraction = 0.0;
        if (regs->is_set(buf)) molefraction = regs->get(buf);

        // Region
        sprintf(buf, "region%d", region_count-1);
        regions.push_back(new Region(buf));
        regions[region_count]->set_index(region_count);
		number_to_materialname(mat_type, buf);
        regions[region_count]->set_material_name(buf);

        logmsg->emit(LOG_INFO,"Region %d (material %s): length=%g, dx=%g...", region_count, buf, length, dx);

        // regions[region_count]->set_material_molefraction(molefraction); LATER (prepare_molefraciton)
        // regions[region_count]->set_material(PropertyContainer<double>* mat); LATER (main.cpp)

        vl = vertex_count-1; // first vertex within layer
        double xl = vertices[vl]->get_coordinate(0);

        while (true) {
            double x_tmp = x;

            x += dx;
            if (x-xl > length) {
                x = xl+length; // we have last vertex of region
            }

            if (fabs(x-x_tmp) < constants::convert_from_SI(units::length, 1e-12)) {
                break; // there is nothing more to do
            }

            //cout << "vertex_count=" << vertex_count << ": x=" << x << endl;
            vertices.push_back(new Vertex(vertex_count, x));
            molefraction_per_vertex.push_back(molefraction);
            vertex_count++;

            edges.push_back(new Edge(edge_count, vertices[vertex_count-1], vertices[vertex_count-2]));
            edge_count++;

            elements.push_back(new Element(elem_count, element_type::interval));
            elem_count++;
            elements[elem_count-1]->add_vertex(vertices[vertex_count-1]);
            elements[elem_count-1]->add_vertex(vertices[vertex_count-2]);
            elements[elem_count-1]->add_edge(edges[edge_count-1]);
            elements[elem_count-1]->set_region(regions[region_count]);
            regions[region_count]->add_element(elements[elem_count-1]);
        }

        region_count++;
    }
    NEGF_ASSERT(region_count>=1, "no region was found!");

    // right contact:
    // - add Region object containing 3 vertices, 2 elements, with spacing = last spacing
    // - add Contact object with only 1 vertex: 3rd-last vertex of structure
    contacts[1]->add_vertex(vertices[vertex_count-1]);
    vertices[vertex_count-1]->set_index_internal(-1);
    double dx_last = vertices[vertex_count-1]->get_coordinate(0) - vertices[vertex_count-2]->get_coordinate(0);
    regions.push_back(new Region("contact_1"));
    regions[region_count]->set_index(region_count);
    regions[region_count]->set_material_name(regions[region_count-1]->get_material_name().c_str());
    for (int vv=0; vv<2; vv++) {

        x += dx_last;

        //cout << "(right contact) vertex_count=" << vertex_count << ": x=" << x << endl;
        vertices.push_back(new Vertex(vertex_count, x));
        molefraction_per_vertex.push_back(molefraction_per_vertex[vertex_count-1]);
        vertex_count++;
        vertices[vertex_count-1]->set_index_internal(-1);

        edges.push_back(new Edge(edge_count, vertices[vertex_count-1], vertices[vertex_count-2]));
        edge_count++;

        elements.push_back(new Element(elem_count, element_type::interval));
        elem_count++;
        elements[elem_count-1]->add_vertex (vertices[vertex_count-1]);
        elements[elem_count-1]->add_vertex (vertices[vertex_count-2]);
        elements[elem_count-1]->add_edge   (edges[edge_count-1]);
        elements[elem_count-1]->set_region (regions[region_count]);
        regions [region_count]->add_element(elements[elem_count-1]);
    }
    region_count++;


    // --------------------------------------
    // add everything to a Geometry object
    // --------------------------------------

    Geometry * geom = new Geometry(vertex_count, edge_count, 0, elem_count);
    for (int ii=0; ii<vertex_count; ii++) geom->add_vertex (vertices[ii]);
    for (int ii=0; ii<edge_count  ; ii++) geom->add_edge   (edges[ii]   );
    for (int ii=0; ii<elem_count  ; ii++) geom->add_element(elements[ii]);
    for (int ii=0; ii<region_count; ii++) geom->add_region (regions[ii] );
    for (int ii=0; ii<2;            ii++) geom->add_contact(contacts[ii] );

    // ------------------
    // prepare and verify
    // ------------------

    logmsg->emit(LOG_INFO_L1,"    Geometry::prepare()...");
    geom->prepare();
    logmsg->emit(LOG_INFO_L1,"    Geometry::verify()...");
    geom->verify();

    // ------------------------------
    // prepare molefractions
    // ------------------------------
    logmsg->emit(LOG_INFO_L2,"    Checking if ternary materials are present...");
    // check if a ternary material is present at all
    bool ternary_material_found = this->check_if_ternary_materials_exist(geom);
    if (!ternary_material_found) {// there is nothing to do
        logmsg->emit(LOG_INFO_L2,"    no ternary materials are present.");
        return geom;
    }

    // vertex-based molefraction field is already stored in molefraction_per_vertex

    // get for each region the value of region-interior points out of the value per vertex
    logmsg->emit(LOG_INFO_L2,"    Reducing vertex-based quantities to region-based quantities...");
    vector<double> x_values = this->get_regionwise_constant_xmole(geom, molefraction_per_vertex);

    // assign molefractions to each region
    logmsg->emit(LOG_INFO_L2,"    Assigning mole fractions...");
    this->assign_xmole_to_regions(geom, x_values); // x_values contains molefraction for each region

    return geom;
}


#ifndef NODFISE

#include "DFISEGrdReader.h"
#include "DFISEGrdWriter.h"
#include "DFISEDatReader.h"
#include "DFISEDatWriter.h"
#include "DFISEPltReader.h"
#include "DFISEPltWriter.h"
using namespace dfise;

#include "InputParserDFISE.cpp"

#endif // NODFISE


/** data of dataset i, value j is stored in values[j*num_datasets+i] */
void InputParser::write_xy(const char* filename, const uint num_values, const uint num_datasets,
							const vector<string> & datanames, const vector<units::UnitType> & unittypes,
							const vector<double> & values) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(num_values>0 && num_datasets>0, "there is nothing to be done.");
	NEGF_ASSERT(filename!=NULL, "encountered null  pointer.");
	NEGF_ASSERT(datanames.size()==num_datasets, "num_datasets does not coincide with number of names.");
	NEGF_ASSERT(values.size()==num_values*num_datasets, "values.size()!=num_values*num_datasets");
	NEGF_ASSERT(unittypes.size()==num_datasets, "sizes of unittype and dataname vectors do not coincide.");
	
	logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s...",filename);
	ofstream fout(filename);
	NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
	fout.precision(12);
	fout.setf(ios::right);
	fout << "# NEGF xy-file\n#";
	for (uint ii=0; ii < num_datasets; ii++) {
		fout << setw(18) << datanames[ii];
		if (ii < num_datasets-1) fout << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < num_values; ii++) {
		for (uint jj=0; jj < num_datasets; jj++) {
			fout << setw(18) << values[ii*num_datasets+jj] / constants::convert_from_SI(unittypes[jj], 1.0);
			if (jj+1 < num_datasets) fout << "\t";
		}
		fout << "\n";
	}
	fout.close();
	logmsg->emit(LOG_INFO_L1, " Done.");
);}


/** helper function for prepare_molefractions(), checks if structure consists of ternary materials 
 *  note: this works only when the molefraction was NOT already appended, i.e. when the material name is still
 *  "AlGaAs" and NOT "AlGaAs0.300" etc...   */
bool InputParser::check_if_ternary_materials_exist(const Geometry * const geom) const
{STACK_TRACE(
	bool result = false;
	for (uint ii = 0; ii < geom->get_num_regions(); ii++) 
	{
		const string & matname = geom->get_region(ii)->get_material_name();
		for (uint jj=0; jj < Constants.ternary_names.size(); jj++) {
			if (Constants.ternary_names[jj]==matname) {
				result = true;
				break;
			}
			if (Constants.ternary_names[jj]==matname.substr(0,Constants.ternary_names[jj].size())) {
				logmsg->emit(LOG_INFO,"Attention: material %s is not treated as an %s ternary material!",
						matname.c_str(), Constants.ternary_names[jj].c_str());
			}
		}
	}
	return result;
);}


/** helper function to get regionwise constant molefractions from a vector which stores the molefractions on each vertex
 *  on the region-interior vertices are considered. */
vector<double> InputParser::get_regionwise_constant_xmole(const Geometry * const geom, const vector<double> & xmole_on_vertices) const
{STACK_TRACE(
	NEGF_ASSERT(xmole_on_vertices.size()==geom->get_num_vertices(), "wrong xmole_vertex vector.");
	vector<double> result;
	result.resize(geom->get_num_regions(), -1.0);
	for (uint ii = 0; ii < geom->get_num_vertices(); ii++) 
	{
		const vector<Region *> & adj_regions = geom->get_regions_near(geom->get_vertex(ii));
		if (adj_regions.size()==1) 
		{
			uint reg = adj_regions[0]->get_index();
			if (result[reg] == -1.0) {
				result[reg] = xmole_on_vertices[ii];
			} else {
				NEGF_ASSERT(result[reg] == xmole_on_vertices[ii], "Nonconstant mole fraction within region.");
			}
		}
	}
	for (uint ii = 0; ii < geom->get_num_regions(); ii++) {
		NEGF_FASSERT(result[ii] != -1.0, "The molefraction for region %s was not found.", 
					geom->get_region(ii)->get_name().c_str() );
	}
	return result;	// copying is OK since we have few regions
);}


/** helper function to assign molefractions to each region */
void InputParser::assign_xmole_to_regions(const Geometry * const geom, const vector<double> & xmole) const
{STACK_TRACE(
	NEGF_ASSERT(xmole.size()==geom->get_num_regions(), "wrong xmole vector.");
	for (uint ii = 0; ii < geom->get_num_regions(); ii++) 
	{
		Region * reg = geom->get_region(ii);
		
		// don't do anything for binary materials
		// NOTE: this also ignores any molefractions defined in these regions!
		if (find(Constants.ternary_names.begin(), Constants.ternary_names.end(), reg->get_material_name())
				 == Constants.ternary_names.end()) {
			continue;
		}
			
		NEGF_ASSERT(xmole[ii]!=-1.0, "the mole fraction for a certain region was not found.");
		logmsg->emit(LOG_INFO_L1,"Setting molefraction to %f for region %s (material %s).",
					xmole[ii],reg->get_name().c_str(), reg->get_material_name().c_str());
		reg->set_material_molefraction(xmole[ii]);
	}
);}


void InputParser::write_xE_matrix(const char* filename, const Matd & matrix, const Geometry * xspace, const vector<double> & energies) const throw (Exception *)	
{STACK_TRACE(
	vector<double> xcoord;
	for (uint ii=0; ii < xspace->get_num_internal_vertices(); ii++) {
		xcoord.push_back(xspace->get_vertex(xspace->get_global_vertex_index(ii))->get_coordinate(0));
	}
	negf::write_xE_matrix(filename, matrix, xcoord, energies); 
);}

void InputParser::write_current_matrix(const char* filename, const Matd & matrix, const Geometry * xspace, const vector<double> & energies) const throw (Exception *)	
{STACK_TRACE(
	vector<double> xcoord;
	for (uint ii=0; ii < xspace->get_num_internal_vertices(); ii++) {
		xcoord.push_back(xspace->get_vertex(xspace->get_global_vertex_index(ii))->get_coordinate(0));
	}
	negf::write_current_matrix(filename, matrix, xcoord, energies); 
);}

void InputParser::write_phi_n_p_Ec_Ev_J(const char* filename, const Geometry * xspace, 
						const vector<double> & pot,  const vector<double> & edens, const vector<double> & hdens, 
						const vector<double> & curr, const vector<double> & Ec,    const vector<double> & Ev)
{STACK_TRACE(
	logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s ...",filename);
	NEGF_ASSERT(pot.size()==xspace->get_num_vertices(), "inconsistent vector size (pot)");
	NEGF_ASSERT(edens.size()==xspace->get_num_vertices(), "inconsistent vector size (edens)");
	NEGF_ASSERT(hdens.size()==xspace->get_num_vertices(), "inconsistent vector size (hdens)");
	NEGF_ASSERT(curr.size()==xspace->get_num_vertices(), "inconsistent vector size (curr)");
	NEGF_ASSERT(Ec.size()==xspace->get_num_vertices(), "inconsistent vector size (Ec)");
	NEGF_ASSERT(Ev.size()==xspace->get_num_vertices(), "inconsistent vector size (Ev)");
	
	ofstream fout(filename);
	NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
	fout.precision(12);
	fout.setf(ios::right);
	fout << "% rows: xcoord, phi, edens, hdens, curr, Ec, Ev\n";
	
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << xspace->get_vertex(ii)->get_coordinate(0) << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << pot[ii] << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << edens[ii] << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << hdens[ii] << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << curr[ii] << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << Ec[ii] << "\t";
	}
	fout << "\n";
	for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
		fout << setw(18) << Ev[ii] << "\t";
	}
	fout << "\n";
	fout.close();
	logmsg->emit(LOG_INFO_L1, " Done.");
);}
