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
#include "Bandedge.h"
using namespace negf;

Bandedge::Bandedge(const Geometry * grid_, const BoxMethod * boxmethod_,
				 const MaterialDatabase * db_, 
				 quantities::PhysicalQuantity electrons_or_holes_, 
				 string verts_or_elems_,
				 double temperature_):
	grid(grid_),
	boxmethod(boxmethod_),
	db(db_),
	electrons_or_holes(electrons_or_holes_),
	verts_or_elems(verts_or_elems_),
	temperature(temperature_),
	corrected(false),
	corrected_eqn(NULL),
	strain_correction(NULL),
	strained(false)
{STACK_TRACE(
	NEGF_ASSERT(grid!=NULL && boxmethod!=NULL && db!=NULL, "encountered null pointer.");
	NEGF_ASSERT(electrons_or_holes==quantities::electron_density 
				|| electrons_or_holes==quantities::hole_density,
				"must specify electrons or holes, not something else.");
	NEGF_ASSERT(verts_or_elems=="vertex" || verts_or_elems=="element",
				"must specify vertex or element, not something else.");
	
	this->its_type = quantities::energy;
		
	// ----------------------------------------------
	// read in values defined on elements from file
	// ----------------------------------------------
	this->elem_values.clear();
	elem_values.resize(grid->get_num_elements(), 0.0);
	for (uint ii = 0; ii < grid->get_num_elements(); ii++)
	{
		const PropertyContainer<double> * mat = grid->get_element(ii)->get_region()->get_material();
		
		// assume storage in eV
		NEGF_FASSERT(mat->is_set("valence_band_edge"), "property \"valence_band_edge\" not set in material %s.", mat->get_name().c_str());
		
		double vbedge  = mat->get("valence_band_edge");
		double cbedge = this->get_cbedge(mat,this->temperature,this->db);
				
		double bandedge = 0.0;
		if (electrons_or_holes==quantities::electron_density) {
			bandedge = cbedge;
		} else if (electrons_or_holes==quantities::hole_density) {
			bandedge = vbedge;
		} else {
			NEGF_EXCEPTION("electrons or holes...");
		}
		elem_values[ii] = constants::convert_from_SI(units::energy, constants::SIec * bandedge);
	}
	
	// -----------------------------------------------
	// set current_variable_values
	// -----------------------------------------------
	if (verts_or_elems=="element")
	{ 
		this->number_of_variables = this->grid->get_num_elements();
		this->current_variable_values.resize(this->number_of_variables, 0.0);
		for (uint ii = 0; ii < grid->get_num_elements(); ii++)
		{
			this->current_variable_values[ii] = elem_values[ii];
		}
	} else if (verts_or_elems=="vertex") 
	{
		this->number_of_variables = this->grid->get_num_vertices();
		this->current_variable_values.resize(this->number_of_variables, 0.0);
		
		const double * const * measure  = boxmethod->get_measure();
		for (uint ii = 0; ii < grid->get_num_vertices(); ii++)
		{
			const vector<Element *> & elems_near_vert = grid->get_elems_near(grid->get_vertex(ii));
			
			// check if grid vertex is at an interface of a quantized with an unquantized region
			bool near_quant   = false;
			bool near_barrier = false;
			
			bool take_barrier_contribs_only = false;
			if (near_quant && near_barrier)
				take_barrier_contribs_only = true;
				
			// check if grid vertex is ***partly*** surrounded by oxide
			// in that case, we will only take the value of the semiconductor material
			const vector<Region *> & adj_regs = grid->get_regions_near(grid->get_vertex(ii));
			bool oxide_found = false;
			bool other_found = false;
			for (uint jj = 0; jj < adj_regs.size(); jj++) {
				bool region_is_oxide = false;
				for (uint kk = 0; kk < Constants.oxide_names.size(); kk++) {
					if (adj_regs[jj]->get_material()->get_name()==Constants.oxide_names[kk]) {
						region_is_oxide = true;
						oxide_found = true;
						break;
					}
				}
				if (!region_is_oxide) other_found = true;
			}
			bool vertex_partly_surrounded_by_oxide = (oxide_found && other_found);
			if (vertex_partly_surrounded_by_oxide) {
				logmsg->emit(LOG_INFO_L2, "vertex %d(%e,%e) is partly surrounded by oxide.",ii,
					grid->get_vertex(ii)->get_coordinate(0),
					grid->get_dimension()>1 ? grid->get_vertex(ii)->get_coordinate(1) : 0.0);
			}
			
			double total_voronoi = 0.0;
			for (uint jj = 0; jj < elems_near_vert.size(); jj++) 
			{
				uint elem_idx = elems_near_vert[jj]->get_index_global();
				uint local_vert_idx = elems_near_vert[jj]->get_local_index(grid->get_vertex(ii));
				
				// special case if grid vertex is *partly* surrounded by oxide
				bool oxide = false;
				if (vertex_partly_surrounded_by_oxide)
				{
					for (uint kk = 0; kk < Constants.oxide_names.size(); kk++) {
						if (elems_near_vert[jj]->get_region()->get_material()->get_name()==Constants.oxide_names[kk]) {
							oxide = true;
						}
					}
				}
				if (oxide) continue;
				
				this->current_variable_values[ii] += measure[elem_idx][local_vert_idx] * elem_values[elem_idx];
				
				total_voronoi += measure[elem_idx][local_vert_idx];
			}
			NEGF_ASSERT(total_voronoi > 0.0, "something went wrong.");
			this->current_variable_values[ii] = current_variable_values[ii] / total_voronoi;
		}
	} else {
		NEGF_EXCEPTION("location must be element or vertex.");
	}
	
	this->timestamp = 999999;	// makes sure eqn is never computed
	
	this->corrected_values.clear();
);}


double Bandedge::get_element_value(uint elem_idx) const
{STACK_TRACE(
	NEGF_ASSERT(elem_values.size() > elem_idx, "invalid element index.");
	return elem_values[elem_idx];
);}


void Bandedge::assign_corrected(Equation * corrected_eqn_)
{STACK_TRACE(
	NEGF_ASSERT(this->corrected, "can only assign a corrected equation when this equation feels it was corrected.");
	this->corrected_eqn = corrected_eqn_; 
);}

Equation * Bandedge::get_corrected_edge() const 
{STACK_TRACE(
	NEGF_ASSERT(corrected_eqn!=NULL, "null pointer encountered."); 
	NEGF_ASSERT(this->corrected, "a corrected equation was found even though the original bandedge was not corrected.");
	return this->corrected_eqn;
);}


/** Workaround to use bowing parameter instead of interpolation of alpha, beta etc. for ternary materials. <BR>
    function is static --> can also be used in other classes */
double Bandedge::get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/, 
							const MaterialDatabase * const database)
{STACK_TRACE(
	int ternary = -1;
	const string & matname = mat->get_name();
	for (uint ii = 0; ii < Constants.ternary_names.size(); ii++) {
		// attention: material name was appended w/ molefraction at this stage
		if (matname.length()>=Constants.ternary_names[ii].length() && 
			matname.substr(0,Constants.ternary_names[ii].length()).compare(Constants.ternary_names[ii]) == 0) {
			ternary = ii;
			break;
		}
	}
	
	 // valence band edge is correct for both binary and ternary materials
	double vbedge = mat->get("valence_band_edge");
	
	// --------------------------
	// case of ternary material
	// --------------------------
	if (ternary!=-1) {
		const TernaryPropertyContainer<double> * matt = database->get_ternary_material(mat->get_name().c_str());
		double molefraction = matt->get_molefraction();
		
		const PropertyContainer<double> * mat0 = matt->get_first_pure_material();
		const PropertyContainer<double> * mat1 = matt->get_second_pure_material();
		NEGF_FASSERT(mat0->is_set("bandgap_0K"), 	"property \"bandgap_0K\" not found in material %s.", 	mat0->get_name().c_str());
		NEGF_FASSERT(mat0->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat0->get_name().c_str());
		NEGF_FASSERT(mat0->is_set("bandgap_beta"), 	"property \"bandgap_beta\" not found in material %s.", 	mat0->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_0K"), 	"property \"bandgap_0K\" not found in material %s.", 	mat1->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat1->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_beta"), 	"property \"bandgap_beta\" not found in material %s.", 	mat1->get_name().c_str());
		
		double bandgap_0K_mat0 	= mat0->get("bandgap_0K");
		double alpha_mat0 		= mat0->get("bandgap_alpha");
		double beta_mat0 		= mat0->get("bandgap_beta");
		double bandgap_0K_mat1  = mat1->get("bandgap_0K");
		double alpha_mat1 		= mat1->get("bandgap_alpha");
		double beta_mat1 		= mat1->get("bandgap_beta");
		
		double bandgap_bowing = 0.0;
		if (matt->is_set("bandgap_bowing")) {
			bandgap_bowing = matt->get("bandgap_bowing");
			// please note that any molefraction dependence of the bowing parameter itself is treated in TernaryPropertyContainer.h
		}
		
		double bandgap_TK_mat0 =bandgap_0K_mat0 - alpha_mat0*T*T / (beta_mat0+T);
		double bandgap_TK_mat1 =bandgap_0K_mat1 - alpha_mat1*T*T / (beta_mat1+T);
		
		double bandgap_TK_alloy = (1.0-molefraction)*bandgap_TK_mat0 + molefraction*bandgap_TK_mat1
								 - molefraction*(1.0-molefraction) * bandgap_bowing;
								 
		return vbedge + bandgap_TK_alloy;
	} else {
	// --------------
	// normal case
	// --------------
		NEGF_FASSERT(mat->is_set("bandgap_0K"), "property \"bandgap_0K\" not found in material %s.", mat->get_name().c_str());
		NEGF_FASSERT(mat->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat->get_name().c_str());
		NEGF_FASSERT(mat->is_set("bandgap_beta"), "property \"bandgap_beta\" not found in material %s.", mat->get_name().c_str());
		
		double bandgap = mat->get("bandgap_0K");
		double alpha   = mat->get("bandgap_alpha");
		double beta    = mat->get("bandgap_beta");
		// assume storage in Kelvin, Tpar=300K, bandgap = Egap @ 0K !!!
		bandgap += /*alpha*300*300 / (beta+300)*/ - alpha*T*T / (beta+T);
		
		return vbedge + bandgap;
	}
);}

/** Add strain correction to current variable values. Any quantum correction then needs to be recomputed. */
void Bandedge::assign_strain_correction(Equation * correction_)
{STACK_TRACE(
	NEGF_ASSERT(!this->is_strained() && this->strain_correction==NULL, "must assign strain corrections only once!");
	NEGF_ASSERT(correction_!=NULL, "null pointer encountered");
	//NEGF_ASSERT(correction_->get_verts_or_elems()==this->verts_or_elems, "both must be defined on the same entity (vertex/element)");
	//NEGF_ASSERT(correction_->get_e_or_h()==this->electrons_or_holes, "CB/VB inconsistency");
	NEGF_ASSERT(correction_->get_type()==quantities::energy, "inconsistent type");
	this->strain_correction = correction_;
		
	// add to current variable values
	for (uint ii = 0; ii < this->get_num_variables(); ii++) {
		this->current_variable_values[ii] += strain_correction->get_value(ii);
	}
	
	// correction in lowdimensional regions must be done again (if already computed)
	this->corrected = false;
	this->corrected_values.clear();
	
	this->strained = true;
);}

