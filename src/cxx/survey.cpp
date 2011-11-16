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


#include "survey.h"

bool survey::fromSurveyFile(char* filename, bool loadSurvey){
	
	//Check if the specified file exists
	ifstream ifile(filename);
	if(! ifile){
		cout << "Survey parameter file " << filename << " does not exists." << endl;
		return false;
	}else{
		ifile.close();
	}
	
	//read parameter file
	paramfile params(filename);
	
	dataFile 	= params.find<string>("survey");
	nbCols   	= params.find<int>("nbcols");
	colPhi   	= params.find<int>("phi");
	colTheta 	= params.find<int>("theta");
	colZ     	= params.find<int>("z");
	colM		= params.find<int>("mass",0);
	string indice 	= params.find<string>("coords");
	if(indice == "G")
		coords = Galactic;
	if(indice == "E")
		coords = Ecliptic;
	rad      	= params.find<bool>("radian");
	sorted   	= params.find<bool>("sorted");
	redshift 	= params.find<bool>("redshift");
	convfc   	= params.find<float>("convfc");
	h 	= params.find<float>("h");
	omg_l	= params.find<float>("omg_l");
	omg_m	= params.find<float>("omg_m");
	omg_b	= params.find<float>("omg_b");
	wa	= params.find<float>("wa");
	w0	= params.find<float>("w0");
	
	if(loadSurvey){
	
	//read actual data
	ifile.open(dataFile.c_str());
	if(! ifile.is_open()){
		cout << "Survey data file " << dataFile << " does not exists." << endl;
		return false;
	}
	
	int64 nlines=0;
	std::string tempString;
	while ( std::getline(ifile, tempString) )
		++nlines;
	cout << nlines << endl;
	
	ifile.clear(); 
	ifile.seekg(0); 
	
	//allocate data array
	if(colM == 0){
		data.alloc(3,nlines);
	}else {
		data.alloc(4,nlines);
	}

	std::cout << "Reading input file" << std::endl;  
	double temp;
	for(int64 i=0; i < nlines ; ++i){
		for(int j=1; j <= nbCols ; ++j){
			ifile >> skipws >> temp;
			if(j == colPhi){
				data(0,i) = temp;
			}else if( j == colTheta){
				data(1,i) = temp;
			}else if( j == colZ){
				data(2,i) = temp;
			}else if(j == colM) {
				data(3,i) = temp;
			}

		}
	}
	
	ifile.close();
	}
	
}

bool survey::toSurveyFile(char* filename){
	std::stringstream paramsFilename;
	std::stringstream dataFilename;
	
	paramsFilename << filename << "_param.txt";
	dataFilename   << filename << ".surv";
	
	//Open 
	ofstream ofile(paramsFilename.str().c_str());
	if(! ofile.is_open()){
		cout << "Survey parameter file " << paramsFilename << " cannot be created." << endl;
		return false;
	}
	
	//writing parameters
	
	ofile << "survey = " << dataFilename << std::endl;
	ofile << "nbcols = " <<    data.size1()   << std::endl;
	ofile << "phi = "    << 1       << std::endl;
	ofile << "theta = "  << 2    << std::endl;
	ofile << "z = "      << 3         << std::endl;
	if (data.size1() == 4)
		ofile << "mass = "   << 4         << std::endl;
	ofile << "coords = " << (coords == Galactic ? "G" : "E") << std::endl;
	ofile << "radian = " << rad          << std::endl;
	ofile << "sorted = " << sorted       << std::endl;
	ofile << "redshift = " << redshift	 << std::endl;
	ofile << "convc = "    << convfc	 << std::endl;
	ofile << "h = " << h << std::endl;
	ofile << "omg_l = " << omg_l << std::endl;
	ofile << "omg_m = " << omg_m << std::endl;
	ofile << "omg_b = " << omg_b << std::endl;
	ofile << "wa = " << wa << std::endl;
	ofile << "w0 = " << w0 << std::endl;
	
	ofile.close();
	
	//write actual data
	ofile.open(dataFilename.str().c_str());
	if(! ofile.is_open()){
		cout << "Survey data file " << dataFilename << " cannot be created." << endl;
		return false;
	}
	
	ofile.setf(std::ios::scientific);
	ofile.precision(7);
	for (int64 i = 0; i < data.size2(); ++i) {
		
		if (data.size1() == 4) {
			ofile <<  data(0,i) << '\t'<< std::showpos << data(1,i) << '\t' << std::noshowpos <<  data(2,i) << '\t' << data(3,i) << std::endl;
		}else {
			ofile <<  data(0,i) << '\t'<< std::showpos << data(1,i) << '\t' << std::noshowpos <<  data(2,i) << std::endl;
		}
		
	}
	
	ofile.close();
	
}


void survey::applySelectionFunction(selectionFunction<double> &function){

	if (data.size1() == 3) {
		planck_rng rng(time(NULL));
		std::list<int64> index;
		list<int64>::iterator it;
		int64 count;
		
		for (int64 i=0; i < data.size2(); ++i) {
			if (rng.rand_uni() <= function.evaluate(data(2,i))/function.evaluate(0.0)) {
				index.push_back(i);
			}
		}
		
		arr2<double> newdata(3,index.size());
		
		count = 0;
		for (it=index.begin() ; it != index.end(); it++) {
			
			newdata(0,count) = data(0,*it);
			newdata(1,count) = data(1,*it);
			newdata(2,count) = data(2,*it);
			count++;
		}
		
		data = newdata;
		
	}else if (data.size1() == 4) {
	  
		
		for (int64 i=0; i < data.size2(); ++i) {
			data(3,i) = function.evaluate(data(2,i))* data(3,i);
		}
	}
}