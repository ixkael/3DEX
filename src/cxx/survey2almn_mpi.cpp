/*
 *  survey2almn.cpp
 *  MRS3D
 *
 *  Created by François Lanusse on 21/06/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */

#include "survey.h"
#include "almn_survey_tools.h"
#include <mpi.h>

extern "C" {
        void __f3dex_fitstools_MOD_almn2fits(char * almnfile, complex16* almn, double* kln,double* cln,int64* nlmax,int64* nmmax,int64* nnmax,char* header,int64* nlheader );
	void __f3dex_utils_MOD_print_almn(complex16* almn,int64* nllim,int64* nmlim,int64* nnlim,int64* nlmax,int64* nmmax,int64* nnmax,char* txt);
	void __f3dex_transforms_MOD_gen_cln(double *cln, double *kln,int64* nnmax,int64* nlmax,double *rmax);
}

int main(int argc, char* argv[]){

  
	MPI::Init(argc,argv);
	
	int NbProc = MPI::COMM_WORLD.Get_size();
	int ProcId = MPI::COMM_WORLD.Get_rank();
	
	if(ProcId == 0){
	std::cout << std::endl;
	std::cout << "+-------------+" << std::endl;
	std::cout << "| survey2almn |" << std::endl;
	std::cout << "+-------------+" << std::endl;
	std::cout << std::endl;
	std::cout << "Computes the Discrete Spherical Fourier-Bessel Transform from an input galaxy survey using the 3DEX algorithm." << std::endl;
	std::cout << std::endl;
	std::cout << "Usage : survey2almn qln_Table_File input_survey output_almn Nside Rmax Lmax Mmax Nmax" << std::endl;
	std::cout << std::endl;
	}
	
	if(argc != 9){
		cout << "Wrong number of arguments" << endl;
		return 0;	
	}	
	
	
	char *qlnTableFile = argv[1];
	char *surveyFile   = argv[2];
	char *almnFile     = argv[3];
	
	int64  nside = atoi(argv[4]);
	double Rmax  = atof(argv[5]);
	
	int64 lmax	 =  atoi(argv[6]);
	int64 mmax	 =  atoi(argv[7]);
	int64 nmax	 =  atoi(argv[8]);
	
	double* qln;
	complex16* almn;
	
	int nmax_local = nmax/NbProc;
	int nOffset = ProcId*nmax_local;
	
	int nmbSurveyCols;
	int64 nmbSurveyObjects;

	survey mySurvey;
	
	if(ProcId == 0){std::cout << "Loading qln...";}
	qln = (double *) malloc((lmax+1)*nmax*sizeof(double));
	if( load_qln(qlnTableFile, nmax, lmax, qln) != 0){
	  std::cout << "The specified qlnFile does not have enough coefficients" << std::endl;
	  return 0;
	}
	if(ProcId == 0){
	  std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Loading survey into memory...";
	}
	//Each machine loads the survey parameters but only proc0 loads the actual survey from the HDD
	for (int i = 0; i < MPI::COMM_WORLD.Get_size(); ++i) {
	    if(MPI::COMM_WORLD.Get_rank() == i){
		if (MPI::COMM_WORLD.Get_rank() == 0) {
			mySurvey.fromSurveyFile(surveyFile);
			nmbSurveyCols = mySurvey.data.size1();
			nmbSurveyObjects = mySurvey.data.size2();
		}else {
		  mySurvey.fromSurveyFile(surveyFile, false);
		}
	    }	
		
		MPI::COMM_WORLD.Barrier();
	}
	if(ProcId == 0){
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	}
	
	// Now Broadcasting size of actual data so that we can allocate survey array for the processes which haven't loaded the survey yet 
	MPI::COMM_WORLD.Bcast(&nmbSurveyCols,1,MPI::INT,0);
	MPI::COMM_WORLD.Bcast(&nmbSurveyObjects,sizeof(int64),MPI::BYTE,0);
	if(MPI::COMM_WORLD.Get_rank() != 0){
	   mySurvey.data.alloc(nmbSurveyCols,nmbSurveyObjects); 
	}
	
	// Now broadcasting the survey data
	if(ProcId == 0){
	std::cout << "Broadcasting survey to every process...";
	}
	MPI::COMM_WORLD.Bcast(&(mySurvey.data(0,0)),sizeof(double)*mySurvey.data.size(),MPI::BYTE,0);
	if(ProcId == 0){
	std::cout << " Done." << std::endl; 
	std::cout << std::endl;
	}
	
	complex16* almnPartiel = (complex16*) malloc((lmax+1)*(mmax + 1)*nmax_local*sizeof(complex16));

	if(ProcId == 0){
	std::cout << "Computing transform..." << std::endl;
	}
	survey2almn(qln, mySurvey, nside,  lmax,  mmax,  nmax_local, Rmax, almnPartiel,nOffset);
	if(ProcId == 0){
	std::cout << " Done." << std::endl; 
	std::cout << std::endl;
	}
	
	if(ProcId == 0){
	    almn = (complex16*) malloc((lmax+1)*(mmax + 1)*nmax*sizeof(complex16));
	    // Copy almns coefficients into output array
#pragma omp parallel for 
		for (int64 l=0 ; l <= lmax ; l++) {
			for (int64 m=0 ; m <= min(mmax,l) ; ++m) {
			   for(int64 n =0; n < nmax_local; n++){
			      almn[m*(lmax +1)*nmax + l*nmax + n] = almnPartiel[m*(lmax +1)*nmax_local + l*nmax_local + n];
			   }
			}
		}
	}
	
	for (int i = 1; i < NbProc; ++i) {
	  if(ProcId == i){
	      MPI::COMM_WORLD.Send(almnPartiel, (lmax+1)*(mmax + 1)*nmax_local*sizeof(complex16), MPI::BYTE, 0, i);
	  }else if(ProcId == 0){
	      MPI::COMM_WORLD.Recv(almnPartiel, (lmax+1)*(mmax + 1)*nmax_local*sizeof(complex16), MPI::BYTE, i, i);
	    
	    // Copy almns coefficients into output array
#pragma omp parallel for 
		for (int64 l=0 ; l <= lmax ; l++) {
			for (int64 m=0 ; m <= min(mmax,l) ; ++m) {
			     for(int64 n =0; n < nmax_local; n++){
			      almn[m*(lmax +1)*nmax + l*nmax + n+ i*nmax_local] = almnPartiel[m*(lmax +1)*nmax_local + l*nmax_local + n];
			   }
			}
		}
	  }
		
	  MPI::COMM_WORLD.Barrier();
	}
	
	
	if (ProcId == 0) {
	  	std::cout << "Saving almn... ";
	int64 nlheader = 120;
	char* header = (char *) malloc(nlheader*80*sizeof(char));
	for(int i=0; i < nlheader*80; i++){
	  header[i] = ' ';
	}
	
	//compute kln
	double *kln = (double *) malloc((lmax+1)*nmax*sizeof(double));
	for(int i=0; i < nmax*(lmax+1); i++){kln[i] = qln[i]/Rmax;}
	
	//compute cln
	double *cln = (double *) malloc((lmax+1)*nmax*sizeof(double));
	__f3dex_transforms_MOD_gen_cln(cln, kln, &nmax, &lmax, &Rmax);
	    
	__f3dex_fitstools_MOD_almn2fits(almnFile, almn, kln, cln, &lmax,&mmax,&nmax, header,&nlheader);
	std::cout << " Done" << std::endl;
	
	std::cout << "Printing the first terms of the decomposition: " << std::endl;
	int64 lim = 15;
	 __f3dex_utils_MOD_print_almn(almn,&lim,&lim,&lim,  &lmax,&mmax,&nmax,"almn");
	 
	 free(almn);
	 free(kln);
	 free(cln);
	}
	
	free(almnPartiel);
	
	MPI::COMM_WORLD.Barrier();
	
	MPI::Finalize();
}
