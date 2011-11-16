/*
 *  sbtools.cpp
 *  FastDSBT
 *
 *  Created by François Lanusse on 06/05/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */

#include "sbtools.h"

#include <fitshandle.h>
#include <iostream>
#include <fstream>

void BesselRoots(int64_t nmax, int64_t lmax, double* qln){
	__f3dex_transforms_MOD_gen_qln(qln,&nmax,&lmax);
}

double sphericalBesselJ(int64_t l, double x){
	double result;
	
	__f3dex_transforms_MOD_bjl(&l, &x, &result);
	
	return result;
}

void write_qln(char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln){
	fitshandle handle;
	handle.create(std::string(tableFileName));
	
	//creation des colonnes 
	std::vector<fitscolumn> cols;
	cols.push_back(fitscolumn(std::string("Spherical Bessel Zeros (n,l) (FORTRAN ARRAY INDEXATION)"),"unknown",1,planckType<double>()));
	
	handle.insert_bintab(cols);
	
	//Parameters
	handle.set_key("NLMAX",(int64) nlmax, "");
	handle.set_key("NNMAX",(int64) nnmax, "");
	
	handle.write_column_raw_void(1, qln, planckType<double>(), (nlmax+1)*nnmax, 0);
	
	handle.close();
}

int load_qln(const char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln){
	
	//Check if the fitsFile exists
	std::ifstream fitsFile;
	fitsFile.open(tableFileName);
	
	if(fitsFile.is_open()){
		fitsFile.close();
		fitshandle handle;
		
		handle.open(std::string(tableFileName));
		
		handle.goto_hdu(2);
		
		if(handle.key_present("NLMAX") && handle.key_present("NNMAX")){
			
			int nlmax_table, nnmax_table; 
			handle.get_key("NLMAX",nlmax_table);
			handle.get_key("NNMAX",nnmax_table);
			
			if(nnmax_table >= nnmax && nlmax_table >= nlmax){
				
				for (int64_t n = 0; n < nnmax; ++n) {
					handle.read_column_raw_void(1, &(qln[(nlmax+1)*n]) ,planckType<double>(), (nlmax+1), (nlmax_table+1)*n);
				}
				
				handle.close();
				return 0;
			}
		}
		
		handle.close();
	}
	
	return -1;
}
