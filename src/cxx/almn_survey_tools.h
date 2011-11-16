/*
 *  almn_survey_tools.h
 *  MRS3D
 *
 *  Created by François Lanusse on 13/04/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */

#ifndef ALMN_SURVEY_TOOLS_H
#define ALMN_SURVEY_TOOLS_H

#include "survey.h"
#include "sbtools.h"
#include <healpix_base.h>
#include <cstring>

#define _USE_MATH_DEFINES
#include <math.h>

struct complex16 {
	double real,imag;
};

extern "C" {
        void __f3dex_transforms_MOD_alnspring2almn(int64* nSide,int64 *nlmax,int64 *nmmax,double* map,complex16* alm1n,double *zbounds);
}


void survey2almn(double *qln, survey &surv, int64 Nside, int64 lmax, int64 mmax, int64 nmax, double Rmax, complex16* almn, int64 nOffset=0){
	double *map;
	double zbounds[2] = {-1,1};
	Healpix_Base base(Nside,RING,SET_NSIDE);
	//memory allocation for fortran code
	complex16 *almnTemp = (complex16*) malloc((lmax+1)*(mmax + 1)*sizeof(complex16));
	
	//double normfactor= sqrt(2.0/M_PI)/Rmax;
	double normfactor=1.0/Rmax;
	
	//memory allocation
	map =(double *) malloc((lmax+1)*base.Npix()*sizeof(double));	
	
	//Creating indexing table
	arr<int64> pixIndex(surv.data.size2());
	
	#pragma omp parallel for
	for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
		pixIndex[i] = base.ang2pix(pointing(surv.data(1,i),surv.data(0,i)));
	}
	
	for (int64 n=0; n < nmax; ++n) {
		int64 actualN = n + nOffset;
		std::cout << "n > " << actualN + 1 << std::endl;
		memset(map,0.0,(lmax+1)*base.Npix()*sizeof(double));
		if (surv.data.size1() == 3) {
		
			#pragma omp parallel for
			for (int64 l=0; l <= lmax; ++l) {
			    for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
				map[pixIndex[i]*(lmax +1) + l] += sphericalBesselJ(l, qln[(lmax+1)*actualN + l] * surv.data(2,i)/Rmax)*normfactor*qln[(lmax+1)*actualN + l];
			    }
			}
		}else {
			#pragma omp parallel for
			for (int64 l=0; l <= lmax; ++l) {
			      for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
					map[pixIndex[i]*(lmax +1) + l] += surv.data(3,i)*sphericalBesselJ(l, qln[(lmax+1)*actualN + l] * surv.data(2,i)/Rmax)*normfactor*qln[(lmax+1)*actualN + l];
				}
			}
		}
				
		__f3dex_transforms_MOD_alnspring2almn(&Nside,&lmax,&mmax,map,almnTemp,zbounds);

		// Copy almns coefficients into output array
#pragma omp parallel for 
		for (int64 l=0 ; l <= lmax ; l++) {
			for (int64 m=0 ; m <= min(mmax,l) ; ++m) {
			    almn[m*(lmax +1)*nmax + l*nmax + n] = almnTemp[m*(lmax +1) + l];
			}
		}
		
 	}	
	
	free(map);
	free(almnTemp);	
}
/*
void survey2almn_mpi(double *qln, survey &surv, int64 Nside, int64 lmax, int64 mmax, int64 nmax_local, double Rmax, complex16* almn, int64 nOffset){
  
  
  	double *map;
	double zbounds[2] = {-1,1};
	Healpix_Base base(Nside,RING,SET_NSIDE);
	//memory allocation for fortran code
	complex16 *almnTemp = (complex16*) malloc((lmax+1)*(mmax + 1)*sizeof(complex16));
	
	//double normfactor= sqrt(2.0/M_PI)/Rmax;
	double normfactor=1.0/Rmax;
	
	//memory allocation
	map =(double *) malloc((lmax+1)*base.Npix()*sizeof(double));	
	
	//Creating indexing table
	arr<int64> pixIndex(surv.data.size2());
	
	#pragma omp parallel for
	for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
		pixIndex[i] = base.ang2pix(pointing(surv.data(1,i),surv.data(0,i)));
	}
	
	for (int64 n=0; n < nmax_local; ++n) {
	  int64 actualN = n + nOffset;
		std::cout << "n > " << n + 1 << std::endl;
		memset(map,0.0,(lmax+1)*base.Npix()*sizeof(double));
		if (surv.data.size1() == 3) {
		
			#pragma omp parallel for
			for (int64 l=0; l <= lmax; ++l) {
			    for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
				map[pixIndex[i]*(lmax +1) + l] += sphericalBesselJ(l, qln[(lmax+1)*actualN + l] * surv.data(2,i)/Rmax)*normfactor*qln[(lmax+1)*n + l];
			    }
			}
		}else {
			#pragma omp parallel for
			for (int64 l=0; l <= lmax; ++l) {
			      for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
					map[pixIndex[i]*(lmax +1) + l] += surv.data(3,i)*sphericalBesselJ(l, qln[(lmax+1)*actualN + l] * surv.data(2,i)/Rmax)*normfactor*qln[(lmax+1)*n + l];
				}
			}
		}
				
		__f3dex_transforms_MOD_alnspring2almn(&Nside,&lmax,&mmax,map,almnTemp,zbounds);

		// Copy almns coefficients into output array
#pragma omp parallel for 
		for (int64 l=0 ; l <= lmax ; l++) {
			for (int64 m=0 ; m <= min(mmax,l) ; ++m) {
			    almn[m*(lmax +1)*nmax + l*nmax + n] = almnTemp[m*(lmax +1) + l];
			}
		}
		
 	}	
	
	free(map);
	free(almnTemp);	
  
  
  
  
  
  
// 	double *map;
// 	int64 lmax;
// 	Healpix_Base base(Nside,RING,SET_NSIDE);
// 	lmax = almn.Lmax();
// 	
// 	double normfactor= sqrt(2.0/M_PI);
// 	
// 	interface.init_kjln(lmax,nmax,Rmax);
// 	
// 	
// 	//memory allocation
// 	map =(double *) malloc((lmax+1)*base.Npix()*sizeof(double));
// 
// 	//Creating indexing table
// 	arr<int64> pixIndex(surv.data.size2());
// 	
// 	#pragma omp parallel for
// 	for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
// 		pixIndex[i] = base.ang2pix(pointing(surv.data(1,i),surv.data(0,i)));
// 	}
// 	
// 
// 	for (int64 n=0; n < almn.Nmax(); ++n) {
// 		int64 actualN = n + almn.Nmax()*nBlockIndex;
// 		std::cout << "n > " << n << std::endl;
// 		memset(map,0.0,(lmax+1)*base.Npix()*sizeof(double));
// 		
// 		if (surv.data.size1() == 3) {
// #pragma omp parallel for
// 			for (int64 l=0; l <= lmax; ++l) {
// 			  for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
// 				map[pixIndex[i]*(lmax +1) + l] += interface.jlnk(l,actualN,surv.data(2,i),Rmax)*normfactor;
// 			}
// 		}
// 	}else {
// #pragma omp parallel for
// 		for (int64 l=0; l <= lmax; ++l) {
// 			for (int64 i = 0; i < (int64) surv.data.size2(); ++i) {
// 				map[pixIndex[i]*(lmax +1) + l] += surv.data(3,i)*interface.jlnk(l,actualN,surv.data(2,i),Rmax)*normfactor;
// 			}
// 		}
// 	}
// 		
// 		interface.alns2almn(Nside, map, n, almn);
// 	}	
// 	
// 	free(map);	
}*/

#endif
