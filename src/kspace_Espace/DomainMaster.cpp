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
#include "DomainMaster.h"
using namespace negf;


DomainMaster::DomainMaster(DomainNode* node) 
: node_counter(0),
  is_frozen(false)
{	
	this->add_node(node);		
}


/** empty domain master constructor */
DomainMaster::DomainMaster()
: node_counter(0),
  is_frozen(false)
{

} 


// define small comperator
bool index_cmp(DomainPoint* a, DomainPoint* b) {
	return a->get_index() < b->get_index();
}

DomainMaster::DomainMaster(const DomainMaster& copy) 
: root_nodes(copy.root_nodes.size()),
  domain_points(0),
  node_counter(copy.node_counter),  
  is_frozen(copy.is_frozen)
{STACK_TRACE(
	// clone root nodes and collect points
	for(unsigned int ii = 0; ii < copy.root_nodes.size(); ii++) {
		root_nodes[ii] = copy.root_nodes[ii]->clone();
		collect_all_indexed_points(*root_nodes[ii], domain_points);	
	}
			
	// and sort the points
	sort(domain_points.begin(), domain_points.end(), index_cmp);

	NEGF_ASSERT(domain_points.size() == copy.domain_points.size(), "domain_points.size() == copy.domain_points.size()");			
);}

/** deep assignment operator */
const DomainMaster& DomainMaster::operator=(const DomainMaster& rhs) 
{STACK_TRACE(
	if(this == &rhs) {
		return *this;	
	}
	// purge old 
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			NEGF_EXCEPTION("somehow my root node is zero!");	
		}	
	}		
	root_nodes.resize(rhs.root_nodes.size());
  	domain_points.resize(0);
  	node_counter = rhs.node_counter;  
  	is_frozen    = rhs.is_frozen;
	// clone root nodes and collect points
	for(unsigned int ii = 0; ii < rhs.root_nodes.size(); ii++) {
		root_nodes[ii] = rhs.root_nodes[ii]->clone();
		collect_all_indexed_points(*root_nodes[ii], domain_points);	
	}			
	// and sort the points
	sort(domain_points.begin(), domain_points.end(), index_cmp);

	NEGF_ASSERT(domain_points.size() == rhs.domain_points.size(), "domain_points.size() == rhs.domain_points.size()");
	return *this;
);}

DomainMaster::~DomainMaster() 
{STACK_TRACE(
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
	//	NEGF_TRACE("master " << this << " deletes: " << root_nodes[ii]);
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			NEGF_EXCEPTION("somehow my root node is zero!");	
		}	
	}
);}

/** reset domain master */
void DomainMaster::clean() 
{STACK_TRACE(
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
	//	NEGF_TRACE("master " << this << " deletes: " << root_nodes[ii]);
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			NEGF_EXCEPTION("somehow my root node is zero!");	
		}	
	}
	domain_points.resize(0); // not my pointers!
	node_counter = 0;
	is_frozen = false;	
);}
	
void DomainMaster::refine() 
{STACK_TRACE(	
	NEGF_ASSERT(!this->frozen(), "can not refine frozen domain!");
	
	int num_points = domain_points.size();
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		root_nodes[ii]->refine();
	}
	update();
	ostringstream sout;	
	sout << "refined master domain from " 
	     << num_points << " to " << domain_points.size();
	logmsg->emit(LOG_INFO_L2, sout.str().c_str());
);}	

void DomainMaster::update() 
{STACK_TRACE(
	NEGF_ASSERT(!this->frozen(), "can not update frozen domain!");
	
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		this->update_domain_points(*root_nodes[ii]);		
	}
);}

void DomainMaster::collapse() 
{STACK_TRACE(
	NEGF_ASSERT(!this->frozen(), "can not collapse frozen domain!");
	
	vector<DomainNode*>::iterator it = root_nodes.begin();
	while(it != root_nodes.end()) {		 
		if(!(*it)->leaf()) {
			// store pointer (we will later delete the node)
			DomainNode* pnode = (*it);
			// remove node (this makes the iterator useless)	
			root_nodes.erase(it);
			// collapse its tree 
			this->collapse_node_tree(*pnode);
			// delete it
			delete pnode;
			// start fresh (collapse node tree possibly added new elements)
			it = root_nodes.begin(); 			
		} else {
			it++;
		}
	}
);}

void DomainMaster::collapse_node_tree(DomainNode& node) 
{STACK_TRACE(
	assert(!node.leaf());
	// process 
	while(node.get_number_of_children() > 0) {		
		unsigned int idx = node.get_number_of_children() - 1;
		// -------------------------------------------
		// if child its not a leaf, descend
		// -------------------------------------------		
		if(!node.get_child(idx).leaf()) {
			this->collapse_node_tree(node.get_child(idx));
			// and delete
			delete node.pop_child();
		} else {
			// -------------------------------------------
			// any leafs will be added to global tree
			// -------------------------------------------
			root_nodes.push_back(node.pop_child());			
		}
	}
);}

/** descends tree and collects pointers of all points which are already indexed */
void DomainMaster::collect_all_indexed_points(DomainNode& node, vector<DomainPoint*>& points) 
{STACK_TRACE(
	if(node.leaf()) {
		if(node.get_point().get_index() != -1) {
			points.push_back(&node.get_point());
		}
	} else {
		for(unsigned int ii = 0; ii < node.get_number_of_children(); ii++) {
			collect_all_indexed_points(node.get_child(ii), points);
		}
	}
);}

void DomainMaster::update_domain_points(DomainNode& node) 
{STACK_TRACE(
	// ------------------------------------------------
	// add leafs points if i don't know them already
	// ------------------------------------------------
	if(node.leaf()) {
		DomainPoint& the_point = node.get_point();
		if(the_point.get_index() == -1) {
			the_point.set_index(node_counter++);
			this->domain_points.push_back(&the_point);
		}
	} else {
		for(unsigned int ii = 0; ii < node.get_number_of_children(); ii++) {
			this->update_domain_points(node.get_child(ii));		
		}
	}
);}


unsigned short DomainMaster::get_dimension() const 
{STACK_TRACE(
	NEGF_ASSERT(this->root_nodes.size() > 0, "empty domain master"); 
	return root_nodes.front()->get_dimension(); 
);}


bool DomainMaster::radial() const
{STACK_TRACE(
	NEGF_ASSERT(this->root_nodes.size() > 0, "empty domain master"); 
	return root_nodes.front()->radial();
);}


const DomainPoint& DomainMaster::get_point(uint idx) const 
{STACK_TRACE(
	NEGF_FASSERT(idx < domain_points.size(), "idx(%d) < domain_points.size()(%d)",idx,domain_points.size()); 
	return *domain_points[idx];
);}

double DomainMaster::get_total_weight() const 
{STACK_TRACE(
	double weight = 0.0;
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		weight += root_nodes[ii]->get_weight();	
	}	
	return weight;
);}

		
DomainMaster::point_iterator DomainMaster::begin() {
	return this->domain_points.begin();	
}

DomainMaster::point_const_iterator DomainMaster::begin() const {
	return this->domain_points.begin();	
}
	 
void DomainMaster::add_node(DomainNode* node) 
{STACK_TRACE(
	NEGF_ASSERT(!this->frozen(), "can not add nodes to frozen domain!");
	
	if(node->get_weight() <= 0.0) {
		logmsg->emit(LOG_WARN,"bad node? %e",node->get_weight());
	}
		
	if(this->root_nodes.size() > 0) {
		NEGF_ASSERT(this->root_nodes.back()->get_dimension() == node->get_dimension(), "a domain master may only have domain of equal dimension");
		NEGF_ASSERT(this->root_nodes.back()->radial() == node->radial(), "you can not mix radial and nonradial nodes!");		
	}
	this->root_nodes.push_back(node);
	update();	
);}

const DomainNode& DomainMaster::get_root_node(unsigned int idx) const 
{STACK_TRACE(
	NEGF_ASSERT(idx < root_nodes.size(), "idx < root_nodes.size()");
	return *root_nodes[idx];	
);}

/** freeze domain -> replace all nodes by frozen ones 
 * 
 * so we can be serialized!
 */
void DomainMaster::freeze() 
{STACK_TRACE(
	NEGF_ASSERT(!this->frozen(), "domain master is already frozen!");
	this->update();
	this->collapse();
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		// make a freeze
		DomainNodeFreeze* frosty = new DomainNodeFreeze(*root_nodes[ii]);
		// delete current node
		delete root_nodes[ii];
		// and replace by frozen one
		root_nodes[ii] = frosty;		
	}
	this->is_frozen = true;
);}

DomainMaster::DomainMaster(istream& in)
: node_counter(0),
  is_frozen(true)
{
	this->read_binary(in);	
}

/** writes frozen nodes to stream. expects domain master to be frozen! */
void DomainMaster::write_binary(ostream& out) const 
{STACK_TRACE(
	NEGF_ASSERT(this->frozen(), "writing does only work on frozen domains!");
	NEGF_ASSERT(node_counter == (signed)domain_points.size(), "node_counter == domain_points.size()"); 
	
		
	int magic = 1818;
	out.write((char*)&magic, sizeof(int));
	
	// ------------------------------------------------
	// just store the number of nodes 
	// ------------------------------------------------
	out.write((char*)&node_counter, sizeof(int));
	
	// ------------------------------------------------
	// and let the nodes to their business
	// ------------------------------------------------
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		DomainNodeFreeze* node = dynamic_cast<DomainNodeFreeze*>(root_nodes[ii]);
		NEGF_ASSERT(node != 0, "could not cast node!");
		node->write_binary(out);	
	}
	out.write((char*)&magic, sizeof(int));			
);}

/** read from stream, expects to be an empty class! */
void DomainMaster::read_binary(istream& in) 
{STACK_TRACE(
	NEGF_ASSERT(root_nodes.size() == 0, "reading binary does only work for empty domains! root_nodes.size() == 0");
	node_counter = 0;
  	is_frozen    = true;
	// read magic
	int magic;
	in.read((char*)&magic, sizeof(int));
	NEGF_ASSERT(magic == 1818, "magic key at domain start is wrong!");
	
	// -----------------------------------------------
	// read number of nodes
	// -----------------------------------------------
	in.read((char*)&node_counter, sizeof(int));
	
	ostringstream sout;	
	sout << "DomainMaster: reading " << node_counter << " nodes from stream";	
	logmsg->emit(LOG_INFO_L2, sout.str().c_str());
	
	// -----------------------------------------------
	// read frozen nodes
	// -----------------------------------------------
	for(int ii = 0; ii < node_counter; ii++) {
		root_nodes.push_back(new DomainNodeFreeze(in));
		collect_all_indexed_points(*root_nodes.back(), domain_points);
	}
	
	// -----------------------------------------------
	// sort points
	// -----------------------------------------------
	sort(domain_points.begin(), domain_points.end(), index_cmp);
	
	NEGF_ASSERT(domain_points.size() == root_nodes.size(), "domain_points.size() == root_nodes.size()"); 

	in.read((char*)&magic, sizeof(int));
	NEGF_ASSERT(magic == 1818, "magic key at domain end is wrong!");			
);}

/** sorter according to coordinates */
bool abs_cmp(DomainPoint* a, DomainPoint* b) {
	return a->get_coord_abs() < b->get_coord_abs();
}

/** reindexes all points according to their absolute coordinate value */
void DomainMaster::reindex_points() 
{STACK_TRACE(
	update();	
	sort(domain_points.begin(), domain_points.end(), abs_cmp);
	// ----------------------------------------------
	// visit all points and reindex
	// ----------------------------------------------
	int index = 0;
	for(unsigned int ii = 0; ii < domain_points.size(); ii++) {
		NEGF_ASSERT(domain_points[ii]->get_index() >= 0, "domain_points[ii]->get_index() >= 0");
		domain_points[ii]->set_index(index++);
	} 
	NEGF_ASSERT(index == node_counter, "index == node_counter");
);}


/** checks if points are indexed equally and are close enought to each other to be accepted as the same grid */
bool DomainMaster::compare_points(const DomainMaster& other) const 
{STACK_TRACE(	
	if(other.get_number_of_points() != this->get_number_of_points()) {
		return false;	
	}	
	unsigned int dimension = this->get_dimension();
	double tmp;
	for(unsigned int ii = 0; ii < domain_points.size(); ii++) {
		if(domain_points[ii]->get_index() != other.domain_points[ii]->get_index()) {
			return false;	
		}	
		tmp = 0.0;
		for(unsigned int cc = 0; cc < dimension; cc++) {
			tmp += negf_math::abs(domain_points[ii]->get_coord(cc) - other.domain_points[ii]->get_coord(cc));			
		}
		if(tmp > 1.0e-10) {
			ostringstream sout;
			sout << "domain points " << ii << " differ by " << tmp;
			logmsg->emit(LOG_INFO_L2, sout.str().c_str());
			return false;	
		}
	}	
	return true;	
);}


/** dump domain content */
void DomainMaster::dump() const 
{STACK_TRACE(
	ostringstream out;
	out << "DomainMaster of " << get_dimension() << "D domain.\n";
	if(radial()) {
		out << "We use the radial approximation.\n";	
	}
	out << "There are " << root_nodes.size() << " root nodes and "
	    << domain_points.size() << " indexed points.\n";
	
	for(unsigned int ii = 0; ii < get_number_of_points(); ii++) {
		out << *domain_points[ii] << "\n"; 	
	}		
	cout << out.str();
);}

