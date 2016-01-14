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
#include "BoxMethod.h"

using namespace negf;

BoxMethod::BoxMethod(const Geometry * const geom_):
	geom(geom_),
	measure(NULL),
	coefficient(NULL),
	node_measure(NULL)
{STACK_TRACE(
	NEGF_ASSERT(geom!=NULL, "null pointer encountered.");
	
	// -----------------------------------------------------
	// allocate space for arrays
	// -----------------------------------------------------
    uint nverts = geom->get_num_vertices();
	uint nelems = geom->get_num_elements();
	node_measure               = new double  [nverts];
	coefficient 	= new double* [nelems];
	measure     	= new double* [nelems];
	for(uint ii=0; ii<nelems; ii++) {
		const int el_nverts = geom->get_element(ii)->get_num_vertices();
		const int el_nedges = geom->get_element(ii)->get_num_edges();
		measure[ii] = new double[el_nverts];
		for(int jj=0; jj<el_nverts; jj++) {
			measure[ii][jj] = 0.0;
		}
		coefficient[ii] = new double[el_nedges];
		for(int jj=0; jj<el_nedges; jj++) {
			coefficient[ii][jj] = 0.0;
		}
	}
	
	// --------------------------------------------------------------
	// compute measure, coefficient
	// --------------------------------------------------------------
	this->compute_1d();
	
	// --------------------------------------------------------------
	// compute node measure
	// --------------------------------------------------------------
	this->compute_node_measures();
	
	// --------------------------------------------------------------
	// compare Voronoi measure and volume calculated from elements
	// --------------------------------------------------------------
	double total_volume1 = 0.0;
	double total_volume2 = 0.0;
	for (uint ii=0; ii<geom->get_num_elements(); ++ii){
		double volume1 = geom->get_element(ii)->get_volume();
		NEGF_ASSERT(volume1 > 0.0, "an element's stored volume was zero.");
		
		double volume2 = 0.0;
		for(uint jj=0; jj < geom->get_element(ii)->get_num_vertices(); jj++) {
			volume2 += measure[ii][jj];
		}
		
		total_volume1 += volume1;
		total_volume2 += volume2;
	}
	logmsg->emit(LOG_INFO,"  total volume from element volumes: %18.12e ", total_volume1);
		
	logmsg->emit(LOG_INFO,"  total volume from voronoi boxes:   %18.12e (Delta = %f %%)", 
					total_volume2, 100.*fabs(total_volume1 - total_volume2)/total_volume1);
	NEGF_ASSERT(fabs(total_volume1 - total_volume2)/total_volume1 < 0.1, 
			"severe discrepancy (>10\%) between total voronoi volume and total element volume detected.");	
);}


BoxMethod::~BoxMethod()
{STACK_TRACE(
	for (uint ii=0; ii<geom->get_num_elements(); ii++) {
		if (measure)     delete [] measure[ii];
		if (coefficient) delete [] coefficient[ii];
	}
	if (measure)     delete [] measure;
	if (coefficient) delete [] coefficient;

	if (node_measure) delete [] node_measure;
);}

/** Get the volume of the Voronoi cell around the Vertex with local index local_vert_idx within Element elem_idx which is contained in that Element. */
double BoxMethod::get_measure(uint elem_idx, uint local_vert_idx) const
{STACK_TRACE(
	NEGF_ASSERT(elem_idx < geom->get_num_elements(), "Invalid element index.");
	NEGF_ASSERT(local_vert_idx < geom->get_element(elem_idx)->get_num_vertices(), "Invalid vertex index.");
	return measure[elem_idx][local_vert_idx];
);}


/** Get the Voronoi surface around the Edge with local index local_edge_idx within Element elem_idx which is contained in that Element. */
double BoxMethod::get_coefficient(uint elem_idx, uint local_edge_idx) const
{STACK_TRACE(
	NEGF_ASSERT(elem_idx < geom->get_num_elements(), "Invalid element index.");
	NEGF_ASSERT(local_edge_idx < geom->get_element(elem_idx)->get_num_edges(), "Invalid edge index.");
	return coefficient[elem_idx][local_edge_idx];
);}


/** Get the total Voronoi volume around a Vertex */
double BoxMethod::get_node_measure(uint vert_idx) const
{STACK_TRACE(
	NEGF_ASSERT(vert_idx < geom->get_num_vertices(), "Invalid vertex index.");
	return node_measure[vert_idx];
);}


double BoxMethod::get_fraction_of_surroundment(const Vertex * vertex, const Region * region) const
{STACK_TRACE(
	
	// get all the surrounding elements of the vertex
	const vector<Element *> & surrounding_elements = this->geom->get_elems_near(vertex);
	
	double total_measure = 0.0;	 // will store the entire voronoi box volume
	double region_measure = 0.0; // will store the volume of the v-box which is inside the given region
	
	for (uint jj = 0; jj < surrounding_elements.size(); jj++)
	{
		Element * element = surrounding_elements[jj];
		Region  * element_region = element->get_region();
		
		// find local index of vertex in the element's vertex list
		uint local_index = element->get_local_index(vertex);

		double measure_contrib = this->measure[element->get_index_global()][local_index];
		total_measure += measure_contrib;
		if (element_region==region) {
			region_measure += measure_contrib;
		}
	}
	
	NEGF_ASSERT(total_measure > 0.0, "something went wrong.");
	return region_measure / total_measure;
);}


double BoxMethod::get_fraction_of_surroundment(const Vertex * vertex, const Element * element) const
{STACK_TRACE(
	
	// get all the surrounding elements of the vertex
	const vector<Element *> & surrounding_elements = this->geom->get_elems_near(vertex);
	
	double total_measure = this->get_node_measure(vertex->get_index_global());
	NEGF_ASSERT(total_measure > 0.0, "something went wrong."); 
	// stores the entire voronoi box volume
	
	vector<Element *>::const_iterator it = find(surrounding_elements.begin(), surrounding_elements.end(), element);
	
	if (it==surrounding_elements.end()) {	// element is not directly adjacent to vertex
		return 0.0;
	} else {
		Element * elem = *it;
		
		// find local index of vertex in the element's vertex list
		uint local_index = elem->get_local_index(vertex);
		
		double element_measure = measure[elem->get_index_global()][local_index]; 
		// stores the volume of the v-box which is inside this element
		
		return element_measure / total_measure;
	}
);}

/** Compute total Voronoi measure and measure within non-oxide materials for each vertex
 *  make sure that this method is called from the derived classes as soon as "measure"-array was set up */
void BoxMethod::compute_node_measures()
{STACK_TRACE(
	for (uint ii=0; ii < geom->get_num_vertices(); ii++) {
		this->node_measure[ii] = 0.0;
	}
	for (uint ii=0; ii<geom->get_num_elements(); ii++) 
	{
		Element * element = geom->get_element(ii);
		
		for (uint jj=0; jj<element->get_num_vertices(); jj++) {
			this->node_measure[element->get_vertex(jj)->get_index_global()] += this->measure[ii][jj];
		}
	}
);}

/** Computes measure and coefficient. Used in 1D and 2D only. */
bool BoxMethod::compute_1d()
{STACK_TRACE(
	NEGF_ASSERT(geom->get_dimension()==1, "expected 1D geometry!");
	
	double dim_scale = 1.0;

	for (uint ii=0; ii<geom->get_num_elements(); ii++)
	{
		Element * element     = geom->get_element(ii);

		for (int jj=0; jj<element->get_num_vertices(); jj++) {
			measure[ii][jj] = 0;
		}
			
		double edge[12]; // maximum for DES_c_Cuboid
		
		NEGF_ASSERT(element->get_type()==element_type::interval, "expected interval.");
		edge[0] = 1;
			
		for (int jj=0; jj<element->get_num_edges(); jj++)
		{
			double length = element->get_edge(jj)->get_length();
			double width  = length * 0.5 * dim_scale * edge[jj];

			// find the local index of the edge's vertices in the element vertex list
			Vertex * v1 = element->get_edge(jj)->get_lower_vertex();
			Vertex * v2 = element->get_edge(jj)->get_upper_vertex();
			int v1_local_index = -1;
			int v2_local_index = -1;
			for (int kk = 0; kk<element->get_num_vertices(); kk++) {
				if (element->get_vertex(kk) == v1) v1_local_index = kk;
				if (element->get_vertex(kk) == v2) v2_local_index = kk;
			}
			NEGF_ASSERT(v1_local_index!=-1 && v2_local_index!=-1, "One of the vertices' local ID was not found");

			// compute measure, coefficient
			this->coefficient[ii][jj]            = edge[jj] / length;
			this->measure[ii][ v1_local_index ] += width;
			this->measure[ii][ v2_local_index ] += width;
		}
	}

	return true;
);}
