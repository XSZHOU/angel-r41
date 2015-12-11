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
#ifndef DOMAINNODE_H_
#define DOMAINNODE_H_

#include "all.h"

#include "DomainPoint.h"

namespace negf {

	/** base class for a collection of vertices (encapsulated, recursive approach) */
	class DomainNode {
		
		friend class DomainNodeFreeze;
		
	public:	
		
		virtual ~DomainNode();
		
		// ---------------------------------------
		// general function
		// ---------------------------------------
		virtual void 		refine() = 0;
		unsigned short 		get_dimension() const { return dimension; }
		virtual DomainNode* clone() const = 0;
		/** returns true if domain is approximated radially */
		virtual bool 		radial() const = 0;
		
		// ---------------------------------------
		// tree functions
		// ---------------------------------------
		uint      			get_number_of_children() const { return children.size(); }
		DomainNode &       	get_child(uint idx) { return const_cast<DomainNode&>(static_cast<const DomainNode&>(*this).get_child(idx)); }
		/** return child */
		const DomainNode & 	get_child(uint idx) const { NEGF_ASSERT(idx < children.size(), "idx < children.size()"); return *children[idx]; }
		/** pop child! returns child with largest index and pops pointer from vector. so you have to delete it later! */
		DomainNode * 		pop_child();
		 	
		// ----------------------------------------
		// leaf functions
		// ----------------------------------------
		/** return true if we are a leaf */
		bool 				leaf() const      { NEGF_ASSERT(point != 0 || children.size() > 0, "point != 0 || children.size() > 0"); return point != 0; }
		DomainPoint & 		get_point()       { NEGF_ASSERT(leaf(), "stupid man, only leaf has a point!"); return *point; }
		/** get leafs point */
		const DomainPoint & get_point() const { NEGF_ASSERT(leaf(), "stupid man, only leaf has a point!"); return *point; }
		/** return weight of domain node (== sum off all child nodes) */
		const double       	get_weight() const { return weight; }
		
		const DomainNode& operator=(const DomainNode& rhs) { NEGF_EXCEPTION("sorry, this is not allowed!"); }
		
	protected:
	
		explicit DomainNode(unsigned short dimension);	// can only be called from derived classes
		explicit DomainNode(const DomainNode& copy);
		vector<DomainNode*>  children;
		DomainPoint*         point;
		double               weight;
		unsigned short       dimension;
		
	};
	
	
	/** frozen arbitrary domain node 
	 *
	 * when a domain master freezes (so, no refinements or updates are possible),
	 * then it collapes and replaces all of its domain nodes by this one
	 *  
	 * the idea is that if i don't refine anymore, i don't need to 
	 * know what i am representing. so i can get rid of the details
	 * 
	 * the advantage: i just have to think about one way of storing
	 * bandstructures and don't have to think how to store the different
	 * types of domains
	 */
	class DomainNodeFreeze : public DomainNode {
	public:
		explicit DomainNodeFreeze(DomainNode& take_over);
		explicit DomainNodeFreeze(istream& in); 
		virtual ~DomainNodeFreeze();
		
		// ---------------------------------------
		// general function
		// ---------------------------------------
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return is_radial; }
		// -----------------------------------------
		// serialization
		// -----------------------------------------
		void write_binary(ostream& out) const;
		
	protected:
		explicit DomainNodeFreeze(const DomainNodeFreeze& copy);
		
	private:
		void read_binary(istream& in);
		bool is_radial;	
	};
	
	
	
	// ====================================================================
	// IMPLEMENTATIONS OF DIFFERENT DIMENSIONS, COORDINATE SYSTEMS
	// ====================================================================
	
	/** line 1D domain (just to be consistent ;-)) */
	class DomainNodeLine : public DomainNode {
	public:
		DomainNodeLine(const double& lx, const double& rx);
		virtual ~DomainNodeLine();	
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return false; }
		
	private:
		DomainNodeLine(DomainPoint* point, const double coords[2]);
		DomainNodeLine(const double coords[2]);
		explicit DomainNodeLine(const DomainNodeLine& copy);
		void init(const double coords[2]);
		/** left x, right x */
		double coords[2];	
	};
	
	/** rectangular 2D domain */
	class DomainNodeRectangle : public DomainNode {
	public:	
		DomainNodeRectangle(const double& llx, const double& lly, const double& urx, const double& ury);
		virtual ~DomainNodeRectangle();	
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return false; }
	private:
		DomainNodeRectangle(DomainPoint* point, const double coords[4]);
		DomainNodeRectangle(const double coords[4]);
		explicit DomainNodeRectangle(const DomainNodeRectangle& copy);
		void init(const double coords[4]);
		/** lower left x, lower left y, upper right x, upper right y */
		double coords[4];
	};
	
	/** Radial approximation of a 2D domain
	 *  an integral over an area, int dx dy, written in the radial approximation is: 2pi * int r dr
	 *  where the integral goes from r0 to r1  
	 */
	class DomainNodeRadialPlane : public DomainNode {
	public:
		DomainNodeRadialPlane(const double& dir_x, const double& dir_y, const double& r0, const double& r1);
		virtual ~DomainNodeRadialPlane();	
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return true; }
	private:
		DomainNodeRadialPlane(DomainPoint* point, const double& dir_x, const double& dir_y, const double& r0, const double& r1);
		explicit DomainNodeRadialPlane(const DomainNodeRadialPlane& copy);
		void init();
		
		/** direction of radial plane */
		double dir_x, dir_y;
		/** interval on radial vector */
		double r0, r1;	
	};
	
	/** cubic 3D domain */
	class DomainNodeCuboid : public DomainNode {
	public:	
		DomainNodeCuboid(const double& llx, const double& lly, const double& llz, const double& urx, const double& ury, const double& urz);
		virtual ~DomainNodeCuboid();	
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return false; }
			
	private:
		DomainNodeCuboid(DomainPoint* point, const double coords[6]);
		DomainNodeCuboid(const double coords[6]);
		explicit DomainNodeCuboid(const DomainNodeCuboid& copy);
		void init(const double coords[6]);
		double coords[6];	
	};
	
	/** singular point with arbitrary weight ... houston, your control! */
	class DomainNodeSingularPoint : public DomainNode {
	public:
		DomainNodeSingularPoint(const double& weight_); // 0D point (used for Quantum Dot Bandstructure)
		DomainNodeSingularPoint(bool radial, const double& weight_, const double& x);
		DomainNodeSingularPoint(bool radial, const double& weight_, const double& x, const double& y);
		DomainNodeSingularPoint(bool radial, const double& weight_, const double& x, const double& y, const double& z);
		virtual ~DomainNodeSingularPoint();	
		virtual void refine();
		virtual DomainNode* clone() const;
		virtual bool radial() const { return is_radial; }
		
	private:
		bool is_radial;
		explicit DomainNodeSingularPoint(const DomainNodeSingularPoint& copy);				
	};

} // namespace negf

#endif /*DOMAINNODE_H_*/
