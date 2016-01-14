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
#include "DomainNode.h"
using namespace negf;


// =========================================================
// DomainNode implementations
// =========================================================


/** protected domain node constructor */
DomainNode::DomainNode(unsigned short dimension_)
: point(0),
  weight(0),
  dimension(dimension_)
{
	//NEGF_TRACE("created DomainNode: " << this);
}

/** copy constructor to perform deep copy */
DomainNode::DomainNode(const DomainNode& copy) 
: point(0),
  weight(copy.weight),
  dimension(copy.dimension)
{
	
	//NEGF_TRACE("copy DomainNode: " << this << " from " << &copy);
	// copy point if necessary
	if(copy.point != 0) {
		point = new DomainPoint(*copy.point);	
	} 
	// copy children
	for(unsigned int ii = 0; ii < copy.children.size(); ii++) {
		children.push_back(copy.children[ii]->clone());
	}
}

/** destructor, deletes children and own point */
DomainNode::~DomainNode() {
//	NEGF_TRACE("deleting DomainNode: " << this);
	for(unsigned int ii = 0; ii < children.size(); ii++) {
		delete children[ii];	
	}
	if(point != 0) {
		delete point; point = 0;	
	}	 	
}

/** pop child from children stack
 * 
 * returns the pointer to the domain node which you then have to delete when 
 * you don't need it anymore
 */
DomainNode* DomainNode::pop_child() {
	if(this->children.size() > 0) {
		DomainNode* tmp = this->children.back();
		this->children.pop_back();
		return tmp; 	
	}	
	return 0;
}




// ===========================================
// implementation of freezed node
// ===========================================

/** take a domain node, steal its point and its weight */
DomainNodeFreeze::DomainNodeFreeze(DomainNode& take_over) 
: DomainNode(take_over.dimension),
  is_radial(false)
{STACK_TRACE(
	NEGF_ASSERT(take_over.leaf(), "Domain Freezing does only work on collapsed trees! every node must be a leaf and the node you passed is not a leaf!");
	NEGF_ASSERT(take_over.point != 0, "damn, the domain node i'm freezing should have a point");
	NEGF_ASSERT(take_over.children.size() == 0, "take_over.children.size() == 0");
	NEGF_ASSERT(negf_math::abs(take_over.weight - take_over.point->get_weight()) < 1.0e-12, "negf_math::abs(take_over.weight - take_over.point->get_weight()) < 1.0e-12");
	this->weight    = take_over.weight;
	this->point     = take_over.point;
	this->is_radial = take_over.radial();
	take_over.point = 0; // steal his point!
);}

/** copy constructor, deep copy */
DomainNodeFreeze::DomainNodeFreeze(const DomainNodeFreeze& copy)
: DomainNode(copy.dimension) 
{
	this->weight    = copy.weight;
	this->point     = new DomainPoint(*copy.point);
	this->is_radial = copy.is_radial;
}

/** initialize node from binary stream */
DomainNodeFreeze::DomainNodeFreeze(istream& in)
: DomainNode(1) // dummy default 
{
	this->read_binary(in);
} 

/** empty destructor */
DomainNodeFreeze::~DomainNodeFreeze() {
	
}

/** no refinement on frozen nodes! */	
void DomainNodeFreeze::refine() 
{STACK_TRACE(
	NEGF_EXCEPTION("no, you can not refine freezed nodes!");
);}

/** clone node */
DomainNode* DomainNodeFreeze::clone() const {
	return new DomainNodeFreeze(*this);	
}

/** write node data to binary stream */
void DomainNodeFreeze::write_binary(ostream& out) const 
{STACK_TRACE(
	NEGF_ASSERT(point != 0, "point != 0");
	NEGF_ASSERT(point->get_index() >= 0, "point->get_index() >= 0");
		
	struct serialdata {
		unsigned int dimension;
		double weight;
		double x, y, z;
		int point_idx;
		bool is_radial;	
	} serial;
	
	serial.dimension = dimension;
	serial.weight    = weight;
	serial.x         = this->point->get_coord(0);
	serial.y         = this->point->get_coord(1);
	serial.z         = this->point->get_coord(2);
	serial.point_idx = this->point->get_index();
	serial.is_radial = is_radial;
	
	// ------------------------------------------------
	// write struct and eigensolutions
	// ------------------------------------------------
	out.write((char*)&serial, sizeof(serialdata));
);}
	
/** reads domain node data and creates a point */	
void DomainNodeFreeze::read_binary(istream& in) 
{STACK_TRACE(	
	// works only for fresh stuff
	NEGF_ASSERT(point == 0, "point == 0"); 
			
	struct serialdata {
		unsigned int dimension;
		double weight;
		double x, y, z;
		int point_idx;
		bool is_radial;	
	} serial;
	
	in.read((char*)&serial, sizeof(serialdata));
	
	dimension = serial.dimension;
	weight    = serial.weight;
    is_radial = serial.is_radial;
	point     = new DomainPoint(serial.x, serial.y, serial.z);
	point->set_index(serial.point_idx);
	point->set_weight(weight);
);}	


// ===========================================
// DomainNodeLine implementation
// ===========================================

/** public constructor for the 1D line domain node */
DomainNodeLine::DomainNodeLine(const double& lx, const double& rx)
: DomainNode(1)
{	
	double coords_[] = {lx, rx};
	point = new DomainPoint((coords_[0] + coords_[1]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeLine::DomainNodeLine(const DomainNodeLine& copy)
: DomainNode(copy) 
{
	for(unsigned short ii = 0; ii < 2; ii++) {
		coords[ii] = copy.coords[ii];
	}		
}

/** creates a deep clone of my node */
DomainNode* DomainNodeLine::clone() const {
	return new DomainNodeLine(*this);	
}

DomainNodeLine::~DomainNodeLine() {
	
}	

/** protected constructor for the line node that gets the domain point if its parent */
DomainNodeLine::DomainNodeLine(DomainPoint* point_, const double coords_[2])
: DomainNode(1) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the 1D line node that creates its own new constructor */
DomainNodeLine::DomainNodeLine(const double coords_[2])
: DomainNode(1) 
{
	point = new DomainPoint((coords_[0] + coords_[1]) / 2.0);
	this->init(coords_);
}	

void DomainNodeLine::init(const double coords_[2]) {	
	for(unsigned short ee = 0; ee < 2; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = negf_math::abs((coords[1] - coords[0]));		     
	point->set_weight(weight); 		
}

/** split the domain node into three new nodes */
void DomainNodeLine::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 3 line elements
		double new_coords[2];
		double seg_x;
		for(unsigned short dd = 0; dd < 3; dd++) {
			seg_x = double(dd) / 3.0;
			// calculate new rectangle boundaries
			new_coords[0] = coords[0] * (1.0 - seg_x)
				          + coords[1] * seg_x; 	
			new_coords[1] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
				          + coords[1] * (seg_x + 1.0 / 3.0);					          
			// ----------------------------------------
			// passing my point
			// ----------------------------------------
			if(dd == 1) {
				children.push_back(
					new DomainNodeLine(
						point,
						new_coords
					)
				);					
				point = 0; 
			} else {
				children.push_back(
					new DomainNodeLine(						
						new_coords
					)
				);
			}		
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}
	
}



// ===========================================================
// DomainNodeRectangle implementation - Rectangular domain
// ===========================================================

/** public constructor for the rectangular domain node */
DomainNodeRectangle::DomainNodeRectangle(const double& llx, const double& lly, const double& urx, const double& ury)
: DomainNode(2)
{	
	double coords_[] = {llx, lly, urx, ury};
	point = new DomainPoint((coords_[0] + coords_[2]) / 2.0, (coords_[1] + coords_[3]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeRectangle::DomainNodeRectangle(const DomainNodeRectangle& copy)
: DomainNode(copy) 
{
	for(unsigned short ii = 0; ii < 4; ii++) {
		coords[ii] = copy.coords[ii];
	}		
}

/** creates a deep clone of my node */
DomainNode* DomainNodeRectangle::clone() const {
	return new DomainNodeRectangle(*this);	
}

DomainNodeRectangle::~DomainNodeRectangle() {
	
}	

/** protected constructor for the rectangular node that gets the domain point if its parent */
DomainNodeRectangle::DomainNodeRectangle(DomainPoint* point_, const double coords_[4])
: DomainNode(2) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the rectangular node that creates its own new constructor */
DomainNodeRectangle::DomainNodeRectangle(const double coords_[4])
: DomainNode(2) 
{
	point = new DomainPoint((coords_[0] + coords_[2]) / 2.0, (coords_[1] + coords_[3]) / 2.0);
	this->init(coords_);
}	

void DomainNodeRectangle::init(const double coords_[4]) {	
	for(unsigned short ee = 0; ee < 4; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = negf_math::abs((coords[2] - coords[0]) * (coords[3] - coords[1]));	     
	point->set_weight(weight); 		
}

/** split the domain node into nine new nodes */
void DomainNodeRectangle::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 9 rectangles
		double new_coords[4];
		double seg_x, seg_y;
		for(unsigned short dd = 0; dd < 3; dd++) {
			for(unsigned short ee = 0; ee < 3; ee++) {
				seg_x = double(dd) / 3.0;
				seg_y = double(ee) / 3.0;
				// calculate new rectangle boundaries
				new_coords[0] = coords[0] * (1.0 - seg_x)
					          + coords[2] * seg_x; 	
				new_coords[1] = coords[1] * (1.0 - seg_y)
					          + coords[3] * seg_y;
				new_coords[2] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
					          + coords[2] * (seg_x + 1.0 / 3.0);					          
				new_coords[3] = coords[1] * (1.0 - seg_y - 1.0 / 3.0)
					          + coords[3] * (seg_y + 1.0 / 3.0);
				// ----------------------------------------
				// passing my point
				// ----------------------------------------
				if(dd == 1 && ee == 1) {
					children.push_back(
						new DomainNodeRectangle(
							point,
							new_coords
						)
					);					
					point = 0; 
				} else {
					children.push_back(
						new DomainNodeRectangle(						
							new_coords
						)
					);
				}	
			}	
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}	
}




// ================================================================================
// DomainNodeRadialPlane implementation - Radial approximation of a 2D domain
// ================================================================================

/** public constructor for the rectangular domain node */
DomainNodeRadialPlane::DomainNodeRadialPlane(const double& dir_x_, const double& dir_y_, const double& r0_, const double& r1_)
: DomainNode(2),
  dir_x(dir_x_),
  dir_y(dir_y_),
  r0(r0_),
  r1(r1_)
{	
	this->init();	
	this->point = new DomainPoint(((r1 + r0) / 2.0) * dir_x, ((r1 + r0) / 2.0) * dir_y);
	point->set_weight(this->weight); 	
}

/** copy constructor for the node */
DomainNodeRadialPlane::DomainNodeRadialPlane(const DomainNodeRadialPlane& copy)
: DomainNode(copy),
  dir_x(copy.dir_x),
  dir_y(copy.dir_y),
  r0(copy.r0),
  r1(copy.r1)
{
	
}

/** creates a deep clone of my node */
DomainNode* DomainNodeRadialPlane::clone() const {
	return new DomainNodeRadialPlane(*this);	
}

DomainNodeRadialPlane::~DomainNodeRadialPlane() {
	
}	

/** protected constructor for the radial node that gets the domain point from its parent */
DomainNodeRadialPlane::DomainNodeRadialPlane(DomainPoint* point_, const double& dir_x_, const double& dir_y_, const double& r0_, const double& r1_)
: DomainNode(2),
  dir_x(dir_x_),
  dir_y(dir_y_),
  r0(r0_),
  r1(r1_)
{
	this->init();
	this->point = point_;	
	point->set_weight(this->weight);
}

/** calculates the radial weight and rescales the direction */
void DomainNodeRadialPlane::init() {
		
	// 2pi ri dr
	this->weight = constants::pi * ((r1 * r1) - (r0 * r0));
	
	double normalize = sqrt(dir_x * dir_x + dir_y * dir_y);
	dir_x /= normalize;
	dir_y /= normalize;
		 	
}

/** split the domain node into three new nodes */
void DomainNodeRadialPlane::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 3 line elements
		double new_coords[2];
		double seg_x;
		for(unsigned short dd = 0; dd < 3; dd++) {
			seg_x = double(dd) / 3.0;
			// calculate new rectangle boundaries
			new_coords[0] = r0 * (1.0 - seg_x)
				          + r1 * seg_x; 	
			new_coords[1] = r0 * (1.0 - seg_x - 1.0 / 3.0)
				          + r1 * (seg_x + 1.0 / 3.0);					          
			// ----------------------------------------
			// passing my point
			// ----------------------------------------
			if(dd == 1) {
				children.push_back(
					new DomainNodeRadialPlane(
						point,
						dir_x,
						dir_y,
						new_coords[0],
						new_coords[1]
					)
				);					
				point = 0; 
			} else {
				children.push_back(
					new DomainNodeRadialPlane(						
						dir_x,
						dir_y,
						new_coords[0],
						new_coords[1]
					)
				);
			}		
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}
}

// ========================================
// DomainNodeSingularPoint implementation 
// ========================================

/** initialize a singular point for a singular point (0D system) */
DomainNodeSingularPoint::DomainNodeSingularPoint(const double& weight_) 
: DomainNode(0),
  is_radial(false)
{
	weight = weight_;
	point  = new DomainPoint(0.0);
	point->set_weight(weight);
}

/** initialize a singular point for a line (probably for outputting reason, as else you could e.g. miss the zero point ...) */ 
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x)
: DomainNode(1),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x);
	point->set_weight(weight);	
}
/** initialize a singular point for a plane (probably for outputting reason, as else you could e.g. miss the zero point ...) */
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x, const double& y) 
: DomainNode(2),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x,y);
	point->set_weight(weight);	
}

/** initialize a singular point for a space (probably for outputting reason, as else you could e.g. miss the zero point ...) */
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x, const double& y, const double& z)
: DomainNode(3),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x,y,z);
	point->set_weight(weight);	
}

/** copy constructor */
DomainNodeSingularPoint::DomainNodeSingularPoint(const DomainNodeSingularPoint& copy)
: DomainNode(copy),
  is_radial(copy.is_radial)
{
}

/** empty destructur */
DomainNodeSingularPoint::~DomainNodeSingularPoint() {}

/** no refinement happens on a singular point */	
void DomainNodeSingularPoint::refine() {}

/** clone node */
DomainNode* DomainNodeSingularPoint::clone() const {
	return new DomainNodeSingularPoint(*this);	
}


// ===========================================
// 3D cuboid implementation
// ===========================================

/** public cuboid constructor */
DomainNodeCuboid::DomainNodeCuboid(const double& llx, const double& lly, const double& llz, const double& urx, const double& ury, const double& urz)
: DomainNode(3)
{	
	double coords_[] = {llx, lly, llz, urx, ury, urz};
	point = new DomainPoint((coords_[0] + coords_[3]) / 2.0, (coords_[1] + coords_[4]) / 2.0, (coords_[2] + coords_[5]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeCuboid::DomainNodeCuboid(const DomainNodeCuboid& copy)
: DomainNode(copy)
{
	for(unsigned short ii = 0; ii < 6; ii++) {
		coords[ii] = copy.coords[ii];
	}
}

/** creates a deep clone of my node */
DomainNode* DomainNodeCuboid::clone() const {
	return new DomainNodeCuboid(*this);	
}

DomainNodeCuboid::~DomainNodeCuboid() {
	
}	

/** protected constructor for the cuboid node that gets the domain point if its parent */
DomainNodeCuboid::DomainNodeCuboid(DomainPoint* point_, const double coords_[6])
: DomainNode(3) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the cuboid node that creates its own new constructor */
DomainNodeCuboid::DomainNodeCuboid(const double coords_[6])
: DomainNode(3) 
{
	point = new DomainPoint((coords_[0] + coords_[3]) / 2.0, (coords_[1] + coords_[4]) / 2.0, (coords_[2] + coords_[5]) / 2.0);
	this->init(coords_);
}	

void DomainNodeCuboid::init(const double coords_[6]) {	
	for(unsigned short ee = 0; ee < 6; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = negf_math::abs((coords[3] - coords[0]) * (coords[4] - coords[1]) * (coords[5] - coords[2]));	     
	point->set_weight(weight); 		
}

/** split the domain node into 27 new nodes */
void DomainNodeCuboid::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 27 cubes
		double new_coords[6];
		double seg_x, seg_y, seg_z;
		for(unsigned short dd = 0; dd < 3; dd++) {
			for(unsigned short ee = 0; ee < 3; ee++) {
				for(unsigned short ff = 0; ff < 3; ff++) {
					seg_x = double(dd) / 3.0;
					seg_y = double(ee) / 3.0;
					seg_z = double(ff) / 3.0;
					// calculate new cube boundaries
					new_coords[0] = coords[0] * (1.0 - seg_x)
						          + coords[3] * seg_x; 	
					new_coords[1] = coords[1] * (1.0 - seg_y)				
						          + coords[4] * seg_y;
					new_coords[2] = coords[2] * (1.0 - seg_z)
						          + coords[5] * seg_z;					          
					new_coords[3] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
						          + coords[3] * (seg_x + 1.0 / 3.0);					          
					new_coords[4] = coords[1] * (1.0 - seg_y - 1.0 / 3.0)
						          + coords[4] * (seg_y + 1.0 / 3.0);
					new_coords[5] = coords[2] * (1.0 - seg_z - 1.0 / 3.0)
						          + coords[5] * (seg_z + 1.0 / 3.0);					          
					// ----------------------------------------
					// passing my point
					// ----------------------------------------
					if(dd == 1 && ee == 1 && ff == 1) {
						children.push_back(
							new DomainNodeCuboid(
								point,
								new_coords
							)
						);					
						point = 0; 
					} else {
						children.push_back(
							new DomainNodeCuboid(						
								new_coords
							)
						);
					}
				}	
			}	
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}	
}

