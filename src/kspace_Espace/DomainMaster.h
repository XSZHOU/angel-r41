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
#ifndef DOMAINMASTER_H_
#define DOMAINMASTER_H_

#include "all.h"
#include "Logger.h"
#include "Exception.h"

#include "DomainNode.h"

namespace negf {


	/** the whole grid over which is integrated */
	class DomainMaster {
	public:
	
		typedef vector<DomainPoint*>::iterator 		 point_iterator;
		typedef vector<DomainPoint*>::const_iterator point_const_iterator;
		
		DomainMaster();
		explicit DomainMaster(DomainNode* node);
		explicit DomainMaster(istream& in); // deserialize
		DomainMaster(const DomainMaster& copy);
		~DomainMaster();
		
		void          			clean(); 
		void         			refine();
		void        			update();
		void        			reindex_points();
		void         			freeze();
		bool        			frozen() const { return is_frozen; }	
		void         			collapse();
		void        			dump() const;
		bool         			compare_points(const DomainMaster& domain) const;
		
		
		double         			get_total_weight() const;
		unsigned short 			get_dimension() const;
		bool           			radial() const;
		uint   					get_number_of_root_nodes() const { return root_nodes.size(); }
		const DomainNode& 		get_root_node(uint idx) const;
			
		point_iterator 			begin();
		point_const_iterator 	begin() const;
		inline point_const_iterator end() const { return domain_points.end(); }
		uint 					get_number_of_points() const { return domain_points.size(); }
		const DomainPoint& 		get_point(uint idx) const;
		const DomainPoint& 		get_first_point() const { return *domain_points.front(); }
		const DomainPoint& 		get_last_point()  const { return *domain_points.back(); } 
			
		void 					add_node(DomainNode* node);
		
		void 					read_binary(istream& in); 		
		void 					write_binary(ostream& out) const;
		
		const DomainMaster& 	operator=(const DomainMaster& rhs);	
		
	private:
		vector<DomainNode*>   	root_nodes;
		vector<DomainPoint*>  	domain_points;
		int  					node_counter;	
		bool 					is_frozen;
		
		void 					update_domain_points(DomainNode& node);
		void 					collapse_node_tree(DomainNode& node);	
		void 					collect_all_indexed_points(DomainNode& node, vector<DomainPoint*>& points);
		
	};

} // namespace negf

#endif /*DOMAINMASTER_H_*/
