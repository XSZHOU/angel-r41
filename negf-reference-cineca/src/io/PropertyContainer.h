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
#ifndef _PROPERTYCONTAINER_H_NEGF
#define _PROPERTYCONTAINER_H_NEGF

#include "all.h"

#include <boost/lexical_cast.hpp>	// $(INC_PATH) must have boost included!

using namespace std;

namespace negf {

	/** general property container class usable for compuational setup, material etc. */
	template<class T>
	class PropertyContainer
	{
	public:
		
		typedef typename map<string, T>::const_iterator 	map_sT_cit;
		typedef typename map<string, T>::iterator 	    	map_sT_it;
	
		PropertyContainer();						// creates noninitialized property container
		PropertyContainer(const char* conffile);	// calls init_from_file
		PropertyContainer(const PropertyContainer&);
		virtual ~PropertyContainer();
	
		// --------------------------------------
		// setup functions
		// --------------------------------------
		void 		 init_from_file(const char* filename)			throw(Exception*);
		void         read_data_from_file(const char* filename)		throw(Exception*);
		void         clear();
		
		void  		 set(const char* key, T value)					throw(Exception*);
		void  		 set(const string& key, T value)				throw(Exception*);	
		void		 set_id(uint id_);
		void		 set_name(const string & name_) { this->name = name_; }
		
		// --------------------------------------
		// access functions
		// --------------------------------------
		virtual bool valid() 							 const throw(Exception*);
		virtual bool valid_value(const char* key,   const T & value) const throw(Exception*); 
		virtual bool valid_value(const string &key, const T & value) const throw(Exception*); 
		virtual bool valid_key	(const char* key)		 const throw(Exception*);
		virtual bool valid_key	(const string& key)		 const throw(Exception*);
		virtual bool is_set		(const char* key)		 const throw(Exception*);
		virtual bool is_set		(const string& key)		 const throw(Exception*);	
		virtual T  	 get		(const char* key)		 const throw(Exception*);
		virtual T 	 get		(const string& key)		 const throw(Exception*);
		uint 		 get_id() 							 const { return this->id; }
		const string & get_name() 						 const { return this->name; }
					
		
	protected:
		/** property definition container */
		struct PropDef {
			string key;
			T error_min, error_max, warn_min, warn_max, default_value;	
			bool mandatory, default_value_avl;		
			bool inside_error_bounds(T value) {
				return value >= error_min && value <= error_max;			
			}
			bool inside_warn_bounds(T value) {
				return value >= warn_min && value <= warn_max;	
			}
			void clear() {
				key.clear(); error_min = error_max = warn_min = warn_max = default_value = (T)0;
				default_value_avl = false;	
			}
		};
		typedef typename map<string, PropDef>::const_iterator 	map_sPD_cit;
		typedef typename map<string, PropDef>::iterator 	    map_sPD_it;
		
		map<string, T> 			value_map;
		map<string, PropDef> 	valid_properties;
		bool   					initialized;
		string 					configuration_file;
		int    					id;
		string					name;
		
		void check_proper_initialization() const;		
		void set_default_values();	
		void strip_whitespaces(string& str);
	};




/** no pointers ... no cleanup */
template<class T>
PropertyContainer<T>::~PropertyContainer() 
{
}

/** create noninitialized property container 
 *
 * property container value files can also point to the corresponding definition
 * file. therefore you can define a uninitialized property container and try to
 * load a value file and get the corresponding definition from that file. 
 * see read_data_from_file. 
 * */
template<class T>
PropertyContainer<T>::PropertyContainer():
	initialized(false),
	name("unknown")
{
}

/** create property container and read configuration from file */
template<class T>
PropertyContainer<T>::PropertyContainer(const char* filename) 
{
	this->init_from_file(filename);	
	this->name = "unknwown";
}

template<class T>
PropertyContainer<T>::PropertyContainer(const PropertyContainer& copy) 
{
	this->value_map.insert(copy.value_map.begin(), copy.value_map.end());
	this->valid_properties.insert(copy.valid_properties.begin(), copy.valid_properties.end());		
	this->initialized = copy.initialized;
	this->name = copy.name;
}
	
template<class T>	
void PropertyContainer<T>::strip_whitespaces(string& str) 
{
	string tmp;
	for(string::iterator it = str.begin(); it != str.end(); it++) {
		if((*it) != ' ' && (*it) != '\t') {
			tmp.push_back(*it);				
		}			
	}
	str = tmp;
}
	
/** read property class setup from file 
 * 
 * instead of hard coding the required values and possible options, they can be
 * easily read from file. the file format must have the following form:
 * # this is a comment. the # must be the first char. 
 * keyname;(mandatory|optional);error_min;error_max;warn_min;warn_max;(default(if optional is set))
 * keyname is the property key
 * (mandatory|optional) denotes whether the specifiy property must be set or could be mandatory
 * error_(min|max) value range where an exception is thrown if such a value is set
 * warn_(min|max)  value range where only a warning message is written to logmsg, but value is accepted
 * default		   if value is optional, then a default value may be set
 * 
 * @param filename property class config filename
 * */
template<class T>
void PropertyContainer<T>::init_from_file(const char* filename) throw(Exception*) 
{STACK_TRACE(
	ifstream fin;
	ostringstream msg;
	string str_in;
	string str_stripped;
	string item;
	string::size_type loc;
	string::size_type last_loc;
	string::iterator it;
	PropDef prop;
	int pcount;
	int line = 1;	
		
	
	// -------------------------------------------------------------
	// open config file
	// -------------------------------------------------------------
	fin.open(filename);
	
	if(!fin) {
		// did not found file in current path. looking for NEGF_MATERIAL_CNFPATH
		char * confpath = getenv("NEGF_MATERIAL_CNFPATH");
		if(confpath != NULL) {			
			string path(confpath);
			string file = path + (path[path.size() - 1] == '/' ? "" : "/") + filename;	
			//ifstream fin2;	// stupid g++4 Red Hat compiler problem necessitates this
			fin.clear();
			fin.open(file.c_str());
			if(!fin) {
				logmsg->emit(LOG_ERROR, "Cannot open or find property definition file %s (1).", file.c_str());
				NEGF_EXCEPTION("Set NEGF_MATERIAL_CNFPATH to point to the configuration directory");
			}
			msg << "reading property container configuration from file " << file;
			logmsg->emit(LOG_INFO_L3, msg.str().c_str());	

//			delete[] confpath;

		} else {
			msg << "Did not find " << filename << " or NEGF_MATERIAL_CNFPATH environment variable.";
			logmsg->emit(LOG_INFO_L1, msg.str().c_str());	
			logmsg->emit(LOG_ERROR, "Did not find %s or NEGF_MATERIAL_CNFPATH.", filename);
			NEGF_EXCEPTION("Set NEGF_MATERIAL_CNFPATH to point to the configuration directory");
		}
	} else {
		msg << "reading property container config from file " << filename;
		logmsg->emit(LOG_INFO_L3, msg.str().c_str());
	}
		
	// ------------------------------------------------------------
	// read config file
	// ------------------------------------------------------------	
	while(getline(fin,str_in)) {
		// ignore comments
		if(str_in.size() > 0 && str_in[0] != '#') {
			str_stripped.clear();
			// strip whitespaces
			this->strip_whitespaces(str_stripped = str_in);
			if(str_stripped.size() > 0) {
				prop.clear();
				// parse line
				pcount   = 0;
				last_loc = 0;
				while((loc = str_stripped.find(";",last_loc)) != string::npos) {	
					item = str_stripped.substr(last_loc, (loc - last_loc));
					try {
						// ------------------------------------------------------------
						// cast into PropDef
						// ------------------------------------------------------------
						switch(pcount) {
							case 0:
								prop.key = item;
								break;
							case 1:
								if(item == "optional") {
									prop.mandatory = false;
								} else if(item == "mandatory") {
									prop.mandatory = true;	
								} else {
									NEGF_FEXCEPTION("in file %s on line %d: expected keyword \"mandatory\" or \"optional\"", filename, line);
								}
								break;
							case 2:
								prop.error_min = boost::lexical_cast<T>(item);
								break;
							case 3:
								prop.error_max = boost::lexical_cast<T>(item);
								// check if min > max
								if(prop.error_min >= prop.error_max) {
									NEGF_FEXCEPTION("in file %s on line %d: error_min >= error_max ", filename, line);
								}
								break;
							case 4:
								prop.warn_min = boost::lexical_cast<T>(item);
								if(prop.warn_min < prop.error_min) {
									NEGF_FEXCEPTION( "in file %s on line %d: warn_min < error_min ", filename, line);
								}
								break;
							case 5:
								prop.warn_max = boost::lexical_cast<T>(item);
								if(prop.warn_min >= prop.warn_max) {
									NEGF_FEXCEPTION( "in file %s on line %d: warn_min >= warn_max ", filename, line);
								}												
								if(prop.warn_max > prop.error_max) {
									NEGF_FEXCEPTION( "in file %s on line %d: warn_max > error_max ", filename, line);
								}							
								break;
							case 6:
								if(prop.mandatory) {
									NEGF_FEXCEPTION("in file %s on line %d: property is mandatory, therefore default value must not be set", filename, line);	
								}
								prop.default_value     = boost::lexical_cast<T>(item);
								prop.default_value_avl = true;
								break;
						}	
					} catch(boost::bad_lexical_cast &) {
						NEGF_FEXCEPTION( "in file %s on line %d: could not parse argument nr. %d: %s", filename, line, pcount + 1, item.c_str());	
					}
					pcount++;	
					last_loc = loc + 1;	
				}  // end while(loc .. )	
				if(pcount < 6) {
					NEGF_FEXCEPTION( "in file %s on line %d: expected at least 5 arguments, found %d", filename, line, pcount - 1);	
				}
				// ------------------------------------------------------------
				// add to valid properties map
				// ------------------------------------------------------------
				if(this->valid_properties.find(prop.key) == this->valid_properties.end()) {		
					this->valid_properties[prop.key] = prop;
				} else {
					NEGF_FEXCEPTION( "in file %s on line %d: property %s already defined", filename, line, prop.key.c_str());
				}
			} // end if(str_stripped.size() > 0)
		} // end if(str[0] != '#')	
		line++;			
	} // end while(getline)	
	this->initialized 		 = true;
	this->configuration_file = filename;
	this->set_default_values();
	fin.close();	
);} 


template<class T>
void PropertyContainer<T>::check_proper_initialization() const {
	if(!this->initialized) {
		NEGF_EXCEPTION("PropertyContainer<T> was not properly initilalized");	
	}	
}


/** read property containers data from given file
 * 
 * the file syntax is simple:
 * # comment (line must start with #)
 * <keyword> = value
 * only keywords that were set in container definiton file are allowed
 * unknown keywords generate exceptions
 * 
 * container objects can also be initialized upon reading a data file.  
 * in this case, the data file must start with a line defining the path to 
 * the definition file
 * # gainconfdef:  <path container definition file>
 * 
 * if the container object already has been initialized and the definition files 
 * differ, the data file is read according to the loaded and not according to 
 * the referenced configuration file. also, a warning will be written to the log
 * 
 */
template<class T>
void PropertyContainer<T>::read_data_from_file(const char* filename) throw(Exception*)
{
	ifstream fin;
	fin.open(filename);
	if(!fin) {
		logmsg->emit(LOG_ERROR, "error while trying to read %s.",filename);
		NEGF_EXCEPTION("can not open data file");
	} 
	string str_in;
	string str_stripped;
	string key;
	T      value;
	string::size_type loc;
	string::iterator it;
	ostringstream msg;
	int line = 1;
			
	msg << "reading data from file " << filename;				
	logmsg->emit(LOG_INFO_L2, msg.str().c_str());
			
	// ------------------------------------------------------------
	// read definition file
	// ------------------------------------------------------------	
	while(getline(fin,str_in)) {
		// --------------------------------------------------------------------------
		// if we're on the first line, we check for any referenced configuration file
		// --------------------------------------------------------------------------
		if(line == 1) {
			//  check if confdef exists
			string confdef_id = "# confdef:";
			if((loc = str_in.find(confdef_id, 0)) != string::npos) {
				string referenced_conffile = str_in.substr(loc + confdef_id.size());
				this->strip_whitespaces(referenced_conffile);	
				// if we are already initialized, we just check if its the same file and emit a warning
				// if its not the same file. if we are not initialized, we try to initialize us
				if(this->initialized) {				
					if(this->configuration_file.compare(referenced_conffile) != 0) {
						std::ostringstream sout;
						sout << "file " << filename << " references property container definition "
							 << referenced_conffile << " while property container is initialized with container definition file: "
							 << this->configuration_file;
						logmsg->emit(LOG_INFO, sout.str().c_str());
					}
				} else {
					this->init_from_file(referenced_conffile.c_str());
				}
			} else if(this->initialized == false) {
				std::ostringstream sout;
				sout << "  unable to read data file " << filename << " due to missing definition";
				NEGF_EXCEPTION(sout.str().c_str());
			}			
		}
		
		// ignore comments
		if(str_in.size() > 0 && str_in[0] != '#') {
			// strip whitespaces
			this->strip_whitespaces(str_stripped = str_in);
			loc = str_stripped.find("#", 0);	// <ss> new: if comment is on same line
			str_stripped = str_stripped.substr(0, loc);
			if(str_stripped.size() > 0) {
				loc = str_stripped.find("=", 0);
				// check if line is valid
				if(loc == string::npos) {
					NEGF_FEXCEPTION("in file %s on line %d: invalid line detected. use <propertykey> = <value>", filename, line);
				}
				key   = str_stripped.substr(0, loc);
				try {
					value = boost::lexical_cast<T>(str_stripped.substr(loc + 1));
				} catch(boost::bad_lexical_cast &) {
					ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": could not parse string \"" << str_stripped.substr(loc + 1) << "\"";
					NEGF_EXCEPTION(sout.str().c_str());
				}						 
				// check if key is valid
				if(!this->valid_key(key)) {
					NEGF_FEXCEPTION( "in file %s on line %d: %s is not a valid key", filename, line, 	key.c_str());
				}
				// check if value is reasonable
				if(!this->valid_properties[key].inside_error_bounds(value)) {
					std::ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": value " << value 
						 << " for key " << key << " is outside error bounds ]" 
						 << this->valid_properties[key].error_min << ","
						 << this->valid_properties[key].error_max << "[";
					NEGF_EXCEPTION(sout.str().c_str());
				}				
				// emit warning if value is out of warning bounds but accept it
				if(!this->valid_properties[key].inside_warn_bounds(value)) {
					std::ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": value " << value 
						 << " for key " << key << " is unreasonable ]" 
						 << this->valid_properties[key].warn_min << ","
						 << this->valid_properties[key].warn_max << "[";
					logmsg->emit(LOG_WARN,sout.str().c_str());
				}						
				// finally, set value	
				this->value_map[key] = value;
			} // end if(str_stripped ..)
		} // end if(str_in.size() ...)
		line++;
	} // end while
	fin.close();
}

/** erases all defined values in value map and resets the all values to the available default values */
template<class T>
void PropertyContainer<T>::clear() {
	this->value_map.clear();
	this->set_default_values();	
}

/** set all default values */
template<class T>
void PropertyContainer<T>::set_default_values() {	
	map_sPD_it it;
	for(it = this->valid_properties.begin(); it != this->valid_properties.end(); it++) {
		if((*it).second.default_value_avl) {
			this->value_map[(*it).first] = (*it).second.default_value;	
		}	
	}	
}

/** returns true if all mandatory fields have a value */
template<class T>
bool PropertyContainer<T>::valid() const throw(Exception *)
{	
	map_sPD_cit it;
	this->check_proper_initialization();	
	// check all properties and check whether a mandatory field is set
	for(it = this->valid_properties.begin(); it != this->valid_properties.end(); it++) {
		if((*it).second.mandatory) {
			if(!this->is_set((*it).first)) {
				return false;	
			}
		}
	}
	return true;
}

/** return true if value is in valid range
 * 
 * throws exception if invalid key is passed
 * 
 * @param key key as defined in PropertyContainer config file
 * @param value value to test
 * @return true if value is inside ]error_min,error_max[, else false
 *  */
template<class T>
bool PropertyContainer<T>::valid_value(const char* key, const T & value) const throw(Exception *)  {
	string tmp(key);
	return this->valid_value(tmp, value);		
}
template<class T>
bool PropertyContainer<T>::valid_value(const string &key, const T & value) const throw(Exception *) {
	this->check_proper_initialization();
	if(!this->valid_key(key)) {
		NEGF_EXCEPTION("invalid key passed to property container");	
	}
	//if(this->valid_properties[key].error_min < value && this->valid_properties[key].error_max > value) {
	map_sPD_cit it = this->valid_properties.find(key);
	if(it != this->valid_properties.end() && (*it).second.error_min < value && (*it).second.error_max > value) {  
		return true;	
	} else {
		return false;	
	}
}

/** check whether key exists in PropertyContainer
 * 
 * does not imply that the value of key is really set. 
 * @param key any property key for any value inside the container
 * @return true, if the property key is defined in the configuration file
 */
template<class T>
bool PropertyContainer<T>::valid_key(const char* key) const throw(Exception *){
	string tmp(key);
	return this->valid_key(tmp);	
}

template<class T>
bool PropertyContainer<T>::valid_key(const string& key) const throw(Exception *){
	this->check_proper_initialization();
    //
    if(this->valid_properties.find(key) == this->valid_properties.end()) {
        //given key exists in the map? not found
		return false;	
	} else {
		return true;	
	}
}

/** check if any value for property key is set
 * 
 * if param is not a valid key, we throw an exception
 * 
 * @param key valid key as defined in the containers configuration file 
 * @return true if key has either a default value (only possible for optional values) or is set
 */
template<class T>
bool PropertyContainer<T>::is_set(const char* key) const throw(Exception*) {
	string tmp(key);
	return this->is_set(tmp);
}


template<class T>
bool PropertyContainer<T>::is_set(const string& key) const throw(Exception*) 
{
	this->check_proper_initialization();
	if(!this->valid_key(key)) {
		NEGF_FEXCEPTION("invalid key %s requested (%s)", key.c_str(), this->get_name().c_str());
	}
	if(this->value_map.find(key) == this->value_map.end()) {
		return false;	
	} else {
		return true;	
	}	
}

/** return value of valid key
 * 
 * @param valid key as defined in container configuration
 * @return value of key if set. if not set, an exception is thrown
 */
template<class T>
T PropertyContainer<T>::get(const char* key) const throw(Exception*) {
	string tmp(key);
	return this->get(tmp);
}
template<class T>
T PropertyContainer<T>::get(const string& key) const throw(Exception*) 
{		
	this->check_proper_initialization();
	
	if(!this->is_set(key)) {
		NEGF_FEXCEPTION("%s: value for \"%s\" is not available", this->name.c_str(), key.c_str());
	}
	map_sT_cit cit = this->value_map.find(key);
	NEGF_ASSERT(cit!=value_map.end(), "key not found.");
	return (*cit).second;
}

/** set value for key
 * 
 * if the passed value is outside of the defined error bounds, an exception will be thrown ...
 * 
 * @param key valid key defined in container configuration. exception thrown if key is invalid
 */
template<class T>
void PropertyContainer<T>::set(const char* key, T value) throw(Exception*) {
	string tmp(key);
    this->set(tmp, value);
}

template<class T>
void PropertyContainer<T>::set(const string& key, T value) throw(Exception*) {
   // std::cout<<"Here we gogo !"<<(this->valid_key(key))<<std::endl;
	this->check_proper_initialization();
    //*{    S.Z.
	if(!this->valid_key(key)) {
		NEGF_FEXCEPTION("invalid key \"%s\"", key.c_str());
	}
	if(!this->valid_properties[key].inside_error_bounds(value)) {
		NEGF_FEXCEPTION("value for \"%s\" is outside error bounds", key.c_str());			
	}	
	if(!this->valid_properties[key].inside_warn_bounds(value)) {			
		logmsg->emit(LOG_WARN, "value for \"%s\" is out of usual bounds", key.c_str()); 			
	}
    //}*/
	this->value_map[key] = value;
}

/** set identifier
 * 
 * used internally by MaterialDatabase (to set the key) and during the assembly procedure 
 * by ProblemDefinition<T>::calculate_element_matrices() to determine cached interaction 
 * matrices. (so interaction matrices are build only once for each material)
 */ 
template<class T>
void PropertyContainer<T>::set_id(unsigned int id_)
{ 
	this->id = id_; 
}


} // end of namespace negf

#endif /*_PROPERTYCONTAINER_H_NEGF*/
