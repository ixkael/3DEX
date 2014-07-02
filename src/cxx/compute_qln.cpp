/*
 *  gen_qln.cpp
 *  MRS3D
 *
 *  Created by François Lanusse on 18/04/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */

#include "sbtools.h"
#include <iostream>
#include <stdlib.h>


int main(int argc, char* argv[]){
	
	std::cout << std::endl;
	std::cout << "+-------------+" << std::endl;
	std::cout << "| compute_qln |" << std::endl;
	std::cout << "+-------------+" << std::endl;
	std::cout << std::endl;
	std::cout << "Tabulates the Bessel roots for further use with survey2almn." << std::endl;
	std::cout << std::endl;
	std::cout << "Usage : gen_qln qln_Table_File Lmax Nmax" << std::endl;
	std::cout << std::endl;
  
	if(argc != 4){
		std::cout << "Wrong number of arguments" << std::endl;
		return(0);
	}
	
	int lmax = atoi(argv[2]);
	int nmax = atoi(argv[3]);
	
	//Allocate qln array
	double *qln = (double *) malloc((lmax +1)*nmax*sizeof(double));
	
	BesselRoots(nmax,lmax,qln);
	
	write_qln(argv[1],nmax,lmax,qln);
	
	std::cout << "Table written" << std::endl;
}
