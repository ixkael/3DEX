/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  François Lanusse <francois.lanusse@supelec.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef SURVEY_H
#define SURVEY_H

#include <arr.h>
#include <paramfile.h>
#include <trafos.h>
#include <planck_rng.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include "selectionFunction.h"

using namespace std;

class survey {

public:
    arr2<double> data;
  
    /*! Creates an empty survey */
    survey() {}
    
    /*! Loads a survey from a parameter file (if loadSurvey=false only the parameters are read not the actual data)*/
    bool fromSurveyFile(char* filename, bool loadSurvey=true);
	
	/*! Applies a radial selection function to the catalog. */
	void applySelectionFunction(selectionFunction<double> &function);
	
	/*! Saves the survey by creating a parameter file along with the actual data file.*/
	bool toSurveyFile(char* filename);
	
    std::string dataFile;
    int			nbCols, colPhi, colTheta, colZ,colM;
    coordsys	coords;
    bool		rad, sorted, redshift;
    double		convfc, h, omg_l, omg_m, omg_b, wa, w0;
    
private:
   
    
};

#endif // SURVEY_H
