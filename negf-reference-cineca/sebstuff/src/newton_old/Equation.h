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
#ifndef EQUATION_H_NEGF
#define EQUATION_H_NEGF

#include "all.h"

using namespace std;

namespace negf {

	/** This is the stem class for an equation which is part of a system solved by Newton iteration. <BR>
	 *  The basic idea is that there is a 1:1 correspondence variable<->equation,
	 *  so each equation has a "native" variable that it solves for (e.g. Poisson->potential). <BR>
	 *  The equation has a list of direct dependencies on other equations (variables).
	 *
	 *  It shall be up to the user which variables are solved in the Newton method.
	 *  This is specified by the boolean solved_in_newton. Variables not solved in the Newton 
	 *  have to be substituted with Newton solution variables. This can only be done when there
	 *  is an explicit equation for the variable (x=f(a,b,..)).
	 *  The class does the substitution automatically and also calculates the derivatives via 
	 *  the chain rule.
	 *
	 *  The Equation class calculates the Newton function F_x of F_x(x)=0 and the corresponding
	 *  lines of the Jacobian matrix (derivatives). Note that communication with other equations 
	 *  is necessary for this.
	 *
	 *  The Equation class also stores current values of its variables, together with a time stamp
	 *  of how recent the values are.
	*/
	class Equation {

		typedef vector<Equation *>::const_iterator eqcit;
		typedef vector<Equation *>::iterator       eqit;

		friend class ExplicitEquation;  
		// needed beause ExplicitEquation::all_direct_derivatives wants to access Equation::all_direct_derivatives
		friend class ExplicitEquationTester;


		public:

			Equation();		//!< Common constructor for all derived classes - does nothing
			virtual ~Equation() {}

			// -------------------------
			// setup functions
			// -------------------------

			/** Construct recursively the list "all_eqns" of all other equations that this equation depends on (including itself) */
			void				handshake() throw(Exception *);

			/** Construct the list of all newton variables that this equation depends on. The result is newton_var_list.
             *  Assumption: all equations are already contained (uniquely) in all_eqns */
			void				set_newton_var_list() throw(Exception *);

			/** Set the global position of the equation in the Newton solver (only relevant if solved_in_newton=true)*/
			void				set_my_offset(uint theoffset) { offset = theoffset; offset_ready = true; }

			/** Renew the current variable values to the time stamp given by tstamp (if possible).
			 *  The function is virtual as it might be overwritten as a hack. <BR>
			 *  The real work is done by compute_value(uint line)
			 *  @param tstamp the desired timestamp. only if the equation's timestamp is smaller will the values be computed. */
			virtual void		compute_values(uint tstamp) throw(Exception *);

			/** Set the current variables to the argument vector.
			 *  The function is virtual as it might be overwritten as a hack */
			virtual void		set_values(const vector<double> & vec, uint tstamp);

			/** Set the reference size */
			void				set_reference_size(double size_) throw(Exception *);

			/** Set the maximum update per variable */
			void				set_maximum_update(double max_update_) throw(Exception *);

			/** Set the time stamp */
			void				set_timestamp(uint tstamp) { timestamp = tstamp; }

			/** Set the name */
			void				set_name(string name_) { this->name = name_; }

			/** Update the current variable values by the vector -update. <BR>
			 *  The function is virtual as it might be overwritten as a hack */
			virtual void		newton_update(const double * update, uint tstamp);

			//!	Calculate and store the number of nonzeros in the jacobian per variable
			void				set_nonzeros();


			// ---------------------------
			// access functions
			// ---------------------------

			/** Get the physical type (density, potential etc.) of the Equation variables */
			quantities::PhysicalQuantity get_type() 		const { return  this->its_type; }

			/** Get the number of variables */
			uint				get_num_variables() 		const { return  this->number_of_variables; }

			/** Get the value of variable idx. The function is virtual as it might be overwritten by derived classes as a hack. */
			virtual double 		get_value(uint idx) 		const throw(Exception *);

			/** Get all the variable values. The function is virtual as it might be overwritten by derived classes as a hack. */
			virtual const vector<double> &	get_values() 	const { return this->current_variable_values; }

			/** Get the reference size of the Equation which will be used when testing for convergence of the Newton iteration. */
			double				get_reference_size()		const { return this->reference_size; }

			/** Get the timestamp, which is used in keeping track of how up-to-date the stored variable values are. */
			uint				get_timestamp() 			const { return this->timestamp; }

			/** Get the name. */
			const string &		get_name()					const { return this->name; }

			/** Check whether the variables of this Equation are Newton solution variables or not (if not, it must be an ExplicitEquation */
			bool				is_solved_in_newton() 		const { return this->solved_in_newton; }

			/** Check whether the recursive list of dependencies on Newton solution Equations has been set up */
			bool				is_newton_var_list_ready() 	const { return this->newton_var_list_ready; }

			/** Check whether the recursive list of all dependencies has been set up */
			bool				is_all_eqns_ready()			const { return this->all_eqns_ready; }

			/** Check whether the offset, i.e. the starting position of this Equation's variables in the Newton solution vector, has been set. */
			bool				is_offset_ready() 			const { return this->offset_ready; }

			/** Get the list of Newton Equation's on which this Equation depends (directly and indirectly) */
			const vector<Equation *> & get_newton_var_list() const throw(Exception *);

			/** Get the list of all Equation's on which this Equation depends (directly and indirectly) */
			const vector<Equation *> & get_all_eqns() 		const throw(Exception *);

			/** Get this Equation's direct dependencies */
			const vector<Equation *> & get_dependencies() 	const { return this->dependencies; }   // use with care!

			/** Get the offset, i.e. the starting position of this Equation's variables in the Newton solution vector */
			uint				get_offset() const throw(Exception *);

			/** Get the stored number of nonzeros w.r.t. a given variable. <BR>
			 *  This is kind of a subset of the information obtained by get_sparsity(), but it is stored and much faster */
			uint				get_nonzeros(uint line) const;

			/** Get the maximum number of nonzeros encountered in a Jacobian line */
			uint				get_max_nonzeros() const;

			/** Check whether the nonzeros information was yet set up */
			bool				nonzeros_was_set() const { return this->nonzeros_set; }

			/** Compute the newton function for a given variable. <BR>
			 *  TO BE IMPLEMENTED BY THE DERIVED CLASS. */
			virtual double		get_newton_function(uint line) const = 0;

			//! Compute a line of the Jacobian matrix \f$ \partial F^i / \partial x_j \f$ belonging to this Equation.
			void				get_newton_derivatives(uint line, uint & nonzeros, uint indices[], double coeffs[]/*, uint tstamp*/) const throw(Exception *);

			/** Compute a list of indices which entries of a line are nonzero in the Jacobian.
			 *  This works almost like get_newton_derivatices(); however, the information about the sparsity is available
			 *  early on and does not rely on any computed values.
			 *  @param line     which variable (0...number_of_variables-1, i.e. without offset)
			 *  @param nonzeros (result) stores the number of nonzeros in the Jacobian at that line.
			 *  @param indices  (result) stores the column indices of the nonzero entries */
			void				get_sparsity(const uint & line, uint & nonzeros, uint indices[]) const throw(Exception *);

			//! Check the timestamps of dependencies; update them if necessary
			void 				check_dependency_timestamps(uint tstamp) throw(Exception *);
		protected:

			// -----------------------------
			// class variables
			// -----------------------------

			// initialized by Equation, ExplicitEquation or ImplicitEquation or by functions of these classes
			vector<double>		current_variable_values;//!< stores the values
			uint				offset; 				//!< global position in Jacobian
			bool				offset_ready;			//!< was the position in the Jacobian set (default: false; unused when not Newton solution variable)
			vector<Equation *>	all_eqns; 				//!< the handshake list storing *all* necessary Equation's to compute this Equation
			bool				all_eqns_ready;			//!< was all_eqns set up (default: false)
			vector<Equation *>	newton_var_list;		//!< a subset of all_eqns containing only the Newton solution variables
			bool				newton_var_list_ready;	//!< was newton_var_list set up (default: false)
			bool				solved_in_newton;		//!< see also ImplicitEquation, ExplicitEquation
			double				reference_size;			//!< gives "big" or "small"... used in the convergence criterion of the Newton solver. default: 1.0
			double 				max_update;				//!< gives maximum update (w/o sign). default: 1e100
			uint				timestamp;				//!< timestamp of the current variable values
			uint				min_dep_timestamp;		//!< stores minimum of the timestamps of the dependencies
			string				name;					//!< the Equation's name, always handy for identification
			vector<uint>		nonzeros_per_line;		//!< stores the number of nonzeros for each line
			uint				max_nonzeros_per_line;	//!< maximum entry in nonzeros_per_line vector
			bool				nonzeros_set;			//!< was the nonzeros information set up (default: false)

			// to be initialized in the constructor of the derived equation classes
			uint 				number_of_variables;	//!< how many variables
			quantities::PhysicalQuantity	its_type;	//!< physical type of the variables (potential, density, ...)
			vector<Equation *>	dependencies; 			//!< other Equation objects that this Equation directly depends on
			// dependencies can NOT be made vector<const Equation *> !!!!

			// --------------------------------------------------------------------------
			// helper functions for public methods, only used in conjunction with these:
			// --------------------------------------------------------------------------

			//! Add idcs2 to idcs1 without duplicates
			void 				add_idcs(uint & nonzeros1, uint idcs1[], uint & nonzeros2, uint idcs2[]) const;
			//! Add idcs2 and cffs2 to idcs1 and cffs1 without duplicates
			void 				add_coeff(uint & nonzeros1, uint idcs1[], double cffs1[], uint & nonzeros2, uint idcs2[], double cffs2[]) const;
			/** Helper for setting up recursively a list with all needed Equaiton objects (directly and indirectly) */
			void                handshake(vector<Equation *> & the_list);

			//! Get the sparsity pattern of a line w.r.t. a certain newton variable
			void	 			get_sparsity(const uint & line, const Equation * nn, uint & nonzeros1, uint icols[]) const;

			/** Get the sparsity pattern w.r.t a direct dependence (eqn must be in the dependency list).
			 *  Implemented by ImplicitEquation and ExplicitEquation classes */
			virtual void		get_dependence_sparsity(const uint & line, const Equation * eqn, 
										uint & nonzeros, uint icol_parts[]) const = 0;

			//! Get the total derivative w.r.t. a newton solution variable. How to do this depends on whether the equation is given explicitly or implicitly!
			void				newton_derivative(const uint & line, const Equation * newton_var, 
										uint & nonzeros, uint indices[], double coeff[]) const;

			/** Get the partial derivative of the Newton function w.r.t. any dependency (eqn must be in the dependency list).
			 *  TO BE IMPLEMENTED BY THE DERIVED CLASS in the case of ImplicitEquation */
			virtual void		get_newton_direct_derivatives(uint line, const Equation * eqn, 
										uint & nonzeros, uint indices[], double coeff[]) const = 0;

			/** Get the partial derivative of the function x = f(..) w.r.t. any dependency (eqn must be in the dependency list).
			 *  Only valid for ExplicitEquation - will throw error when used with an ImplicitEquation */
			virtual void		all_direct_derivatives(const uint & line, const Equation * newton_var, 
										uint & nonzeros, uint indices[], double coeff[]) const = 0;

			/** Does the work of compute_values(timestamp).
			 *  TO BE IMPLEMENTED BY THE DERIVED CLASS in the case of ExplicitEquation;
             *  Invalid method for ImplicitEquation (computation not possible). */
			virtual double		compute_value(uint line) const = 0;
			
	};

}

#endif /*EQUATION_H_NEGF*/
