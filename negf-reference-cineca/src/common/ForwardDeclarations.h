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
#ifndef FORWARDDECLARATIONS_H_
#define FORWARDDECLARATIONS_H_

/** Needed for SWIG: conversion C++ code --> Python modules
 *  SWIG needs to know the namespaces of the classes 
 *  include this .h-file before all other headers in your .i-file */

namespace negf {
	
	class Logger;
	class Timer;
	class Filenames;
	class Interrupt;
	class Matc;
	class BMatc;
	class Matd;
	class Vecd;
	class MPI;
	
	class Vertex;
	class Edge;
	class Element;
	class Region;
	class Face;
	class Iface;
	class Contact;
	class Geometry;	
	class BoxMethod;
		
	class InputParser;	
	class OutputData;	
	class MaterialDatabase;	
	template<class T> class PropertyContainer;
	template<class T> class TernaryPropertyContainer;
	class Options;
		
	class Equation;	
	class ExplicitEquation;	
	class ImplicitEquation;
	class Poisson;
	class VertexData;	
	class EpsilonRegionwise;
		
	class DomainPoint;
	class DomainNode;
	class DomainMaster;
	class Energies;
	class Kspace;
	
	template<class T> class NEGFObject;
	class SelfEnergy;
	class SelfEnergies;
	class GreenFunction;
	class GreenFunctions;
	class InnerLoop;
	
} 

#endif /*FORWARDDECLARATIONS_H_*/
