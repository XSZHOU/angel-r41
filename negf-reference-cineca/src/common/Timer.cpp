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
#include "Timer.h"

using namespace negf;

Timer::Timer()
{STACK_TRACE(
	clicks.clear();
	clicks.push_back(this->get_time());
	formatted_times.clear();
	formatted_times.push_back(this->get_formatted_timestring());
	clicknames.clear();
	string mystart = "start";
	clicknames.push_back(mystart);
);}

double Timer::get_time()
{STACK_TRACE(
	#ifdef WIN32
		struct _timeb timebuffer;
		_ftime( &timebuffer );
		return timebuffer.time + 1e-3*timebuffer.millitm;
	#else
		timeval tval;
		gettimeofday(&tval, NULL);  
		// tval.tv_sec = seconds since 1.1.1970,   tval.tv_usec = microseconds
		return tval.tv_sec + 1e-6 * tval.tv_usec;
	#endif
);}


Timer::~Timer()
{STACK_TRACE(
);}


void Timer::click(string name)
{STACK_TRACE(
	NEGF_ASSERT(name.length()>0, "name of timer click must not be empty.");
	NEGF_ASSERT(find(clicknames.begin(), clicknames.end(), name) == clicknames.end(),
				"clickname must be unique.");
	clicks.push_back(this->get_time());
	formatted_times.push_back(this->get_formatted_timestring());
	clicknames.push_back(name);
);}


double Timer::get_click(uint ii)
{STACK_TRACE(
	NEGF_ASSERT(ii<clicks.size(), "invalid click index.");
	return clicks[ii];
);}


string Timer::get_clickname(uint ii)
{STACK_TRACE(
	NEGF_ASSERT(ii<clicknames.size(), "invalid click index.");
	return clicknames[ii];
);}
	

uint Timer::get_index(string name)
{STACK_TRACE(
	uint idx = 88888888;
	for (uint ii=0; ii<clicknames.size(); ii++) {
		if (clicknames[ii]==name) {
			idx = ii;
			break;
		}
	}
	NEGF_ASSERT(idx!=88888888, "click name not found.");
	return idx;
);}
	
	
double Timer::get_click(string name)
{STACK_TRACE(
	NEGF_ASSERT(clicks.size()==clicknames.size(), "something is wrong.");
	return clicks[get_index(name)];
);}


double Timer::get_seconds_since_start(string name)
{STACK_TRACE(
	double start = clicks[0];
	double now = clicks[get_index(name)];
	return now - start;
);}


double Timer::get_seconds_since_last_click(string clickname)
{STACK_TRACE(
	uint ii = get_index(clickname);
	NEGF_ASSERT(ii>0, "the specified clickname was the first click.");
	return clicks[ii] - clicks[ii-1];
);}


double Timer::get_seconds_between(string name1, string name2)
{STACK_TRACE(
	double click1 = clicks[get_index(name1)];
	double click2 = clicks[get_index(name2)];
	return click1 - click2;
);}


string Timer::get_formatted_timestring()
{STACK_TRACE(
	#ifdef WIN32
		NEGF_EXCEPTION("Not implemented.");
		return "";
	#else
		timeval tval;
		gettimeofday(&tval, NULL);  
		
		char buf[1000];
		sprintf(buf,"%d%02d%02d_%02d%02d%02d",
					(gmtime(&tval.tv_sec))->tm_year+1900,
					(gmtime(&tval.tv_sec))->tm_mon+1,
					(gmtime(&tval.tv_sec))->tm_mday,
					(gmtime(&tval.tv_sec))->tm_hour,
					(gmtime(&tval.tv_sec))->tm_min,
					(gmtime(&tval.tv_sec))->tm_sec);
		string result(buf);
		return result;
	#endif
);}


string Timer::get_formatted_time(string clickname)
{STACK_TRACE(
	uint ii = get_index(clickname);
	NEGF_ASSERT(formatted_times.size() > ii, "something went wrong.");
	return formatted_times[ii];
);}
