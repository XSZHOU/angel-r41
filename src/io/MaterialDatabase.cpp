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
#include "MaterialDatabase.h"
using namespace std;


namespace negf
{

/** main constructor adds current directory to search path per default */
MaterialDatabase::MaterialDatabase()  throw (Exception *)
{STACK_TRACE(
	paths.push_back(".");
);}


MaterialDatabase::~MaterialDatabase()  throw (Exception *)
{STACK_TRACE(
	for(material_iterator it = materials.begin(); it != materials.end(); it++) {
		if((*it).second != 0) {
			delete (*it).second;
			(*it).second = 0;	
		}	
	}		
	// do nothing w/ ternary_materials since all of those are also contained in materials
);}


/** add path to search path for search material files **/ 	
void MaterialDatabase::add_search_path(const char* path)  throw (Exception *)
{STACK_TRACE(
	string tmp(path);
	if(tmp.size() == 0) {
		NEGF_EXCEPTION("tried to add empty path to search path");			
	}
	paths.push_back(tmp);
);}

	
/** get value for some materials property
 * 
 * @param mat material identifier. 
 * @param key property key as defined in a property container config file
 * @return the property value. if either material parameter can not be found or key is invalid, an exception is thrown
 */ 	
double MaterialDatabase::get(const char* mat,   const char* key) throw (Exception *)
{STACK_TRACE(
	if(this->property_is_set(mat, key))	{
		return this->materials[mat]->get(key);
	}	
	NEGF_FEXCEPTION("material property %s for %s is not set", key, mat);		
);}


/** check whether material properties for key are available 
 * @param key material key as referenced in grid file
 * */
bool MaterialDatabase::material_exists(const char* key) const throw (Exception *)
{STACK_TRACE(
	string tmp(key);
	material_citerator it;
	if((it = this->materials.find(tmp)) != this->materials.end()) {
		return true;	
	} else {
		return false;
	}
);}


/** check whether material properties for key are available 
 * @param key material key as referenced in grid file
 * */
bool MaterialDatabase::ternary_material_exists(const char* key) const throw (Exception *)
{STACK_TRACE(
	string tmp(key);
	ternary_material_citerator it;
	if((it = this->ternary_materials.find(tmp)) != this->ternary_materials.end()) {
		return true;	
	} else {
		return false;
	}
);}


/** check whether material has all mandatory values set and all values are inside given bounds 
 * @param key the material key as defined in grid file
 */
bool MaterialDatabase::material_is_valid(const char* key)  throw (Exception *)
{STACK_TRACE(
	if(!this->material_exists(key)) {
		NEGF_FEXCEPTION("unloaded material %s requested", key);
		return false;	
	}
	return this->materials[key]->valid();	
);}


/** check if a property of a material is set
 * @param mat material key
 * @param key property key
 * @return true if property is set, false if not set, exception is thrown if mat is not available or 
 *         key is not a valid key (as defined in the configuration file)
 */
bool MaterialDatabase::property_is_set(const char* mat, const char* key)  throw (Exception *)
{STACK_TRACE(
	if(this->material_exists(mat)) {
		return this->get_material(mat)->is_set(key);
	} else {
		NEGF_FEXCEPTION("request for %s of material %s is invalid as material was never loaded", key, mat);
	}		
);}


/** return pointer to the property container having the material properties */
PropertyContainer<double>* MaterialDatabase::get_material(const char* key) const  throw (Exception *)
{STACK_TRACE(
	string key_str(key);
	material_citerator it;
	if((it = this->materials.find(key_str)) != this->materials.end()) {
		return (*it).second;
	} else {
		string error_msg = "invalid material requested: "; error_msg.append(key);
		error_msg.append("\npresent materials:");
		for (it = this->materials.begin(); it != this->materials.end(); it++) {
			error_msg.append("\n   ");error_msg.append(it->first);
		}
		NEGF_EXCEPTION(error_msg.c_str());
		return 0;
	}		
);}


PropertyContainer<double>* MaterialDatabase::get_material(const int &id) const  throw (Exception *)
{STACK_TRACE(
	material_citerator it;
	for(it = this->materials.begin(); it != this->materials.end(); it++) {
		if((signed)(*it).second->get_id() == id) {
			return (*it).second;	
		}	
	}	
	NEGF_EXCEPTION("invalid material requested");
	return 0;
);}


/** return pointer to the ternary property container */
TernaryPropertyContainer<double>* MaterialDatabase::get_ternary_material(const char* key) const  throw (Exception *)
{STACK_TRACE(
	string key_str(key);
	ternary_material_citerator it;
	if((it = this->ternary_materials.find(key_str)) != this->ternary_materials.end()) {
		return (*it).second;
	} else {
		string error_msg = "invalid ternary material requested: "; error_msg.append(key);
		error_msg.append("\npresent ternary materials:");
		for (it = this->ternary_materials.begin(); it != this->ternary_materials.end(); it++) {
			error_msg.append("\n   ");error_msg.append(it->first);
		}
		NEGF_EXCEPTION(error_msg.c_str());
		return 0;
	}		
);}


/** tries to load the material from filesystem */
void MaterialDatabase::load_material(const char* name)  throw (Exception *)
{STACK_TRACE(
	list<string>::const_iterator it;
	string file;
	ifstream fin;
	for(it = paths.begin(); it != paths.end(); it++) {
		file = (*it) + ((*it)[(*it).size() - 1] == '/' ? "":"/") + name + ".mat";
		logmsg->emit(LOG_INFO_L3, "Trying to open \"%s\" ...",file.c_str());
		fin.open(file.c_str());
		if(fin) {
			//logmsg->emit(LOG_INFO, "Reading \"%s\" params from %s", name, file.c_str());
			PropertyContainer<double>* prop = new PropertyContainer<double>();
			prop->read_data_from_file(file.c_str());
			prop->set_name(name);
			this->add_material(name, prop);
			fin.close();
			return;				
		} else {
			logmsg->emit(LOG_INFO_L3, "%s not found.", file.c_str());
		}
		fin.close();
		fin.clear();	// important for Red Hat gcc compiler versions!
	}
	ostringstream sout;
	sout << "unable to locate material file " << name << ".mat in search path:\n";
	for(it = paths.begin(); it != paths.end(); it++) {
		sout << (*it) << "\n";	
	}
	NEGF_EXCEPTION(sout.str().c_str());		
);}
//}


/** tries to load a ternary material with a certain molefraction filesystem 
 *  assumption: molefraction has been incorporated into last 5 chars of the name(?)*/
void MaterialDatabase::load_ternary_material(const char* name, const double molefraction)  throw (Exception *)
{//STACK_TRACE(
	list<string>::const_iterator it;
	string file;
	ifstream fin;
	for(it = paths.begin(); it != paths.end(); it++) 
	{
		file = (*it) + ((*it)[(*it).size() - 1] == '/' ? "":"/") + name + ".mat";
		logmsg->emit(LOG_INFO_L3, "Trying to open \"%s\" ...",file.c_str());
		fin.open(file.c_str());
		if(fin) {
			logmsg->emit(LOG_INFO, "Reading \"%s\" params from %s", name, file.c_str());
			TernaryPropertyContainer<double>* prop = new TernaryPropertyContainer<double>(file.c_str(), molefraction);
			ostringstream new_name;
			new_name << name << molefraction;
			prop->set_name(new_name.str().c_str());
			this->add_ternary_material(new_name.str().c_str(), prop);
			//char new_name[1000];
			//sprintf(new_name,"%s%5.3f", name, molefraction);
			//prop->set_name(new_name);
			//this->add_ternary_material(new_name, prop);
			fin.close();
			return;				
		} else {
			logmsg->emit(LOG_INFO_L3, "%s not found.", file.c_str());
		}
		fin.close();
		fin.clear();	// important for Red Hat gcc compiler versions!
	}
	ostringstream sout;
	sout << "unable to locate material file " << name << ".mat in search path:\n";
	for(it = paths.begin(); it != paths.end(); it++) {
		sout << (*it) << "\n";	
	}
	NEGF_EXCEPTION(sout.str().c_str());		
}//);}


void MaterialDatabase::add_material(const string& name, PropertyContainer<double>* prop)  throw (Exception *)
{STACK_TRACE(
	if(this->material_exists(name.c_str())) {
		NEGF_FEXCEPTION("material %s already defined", name.c_str());
		return;	
	}
	prop->set_id((uint) this->materials.size());
	this->materials[name] = prop;
	return;
);}


/** same as add_material, but property container is also added to ternary_materials */
void MaterialDatabase::add_ternary_material(const string& name, TernaryPropertyContainer<double>* prop)  throw (Exception *)
{STACK_TRACE(
	if(this->material_exists(name.c_str())) {
		NEGF_FEXCEPTION("material %s already defined", name.c_str());
		return;	
	}
	prop->set_id((uint) this->materials.size());
	this->materials[name] = prop;
	this->ternary_materials[name] = prop;
	return;
);}

} // end namepsace
