/*
 *  sbtools.h
 *  FastDSBT
 *
 *  Created by François Lanusse on 04/05/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */

#ifndef ST_TOOLS_H
#define ST_TOOLS_H

#include <stdint.h>

extern "C" {
	void __f3dex_transforms_MOD_gen_qln(double* qln, int64_t *nnmax, int64_t *nlmax);
	void __f3dex_transforms_MOD_bjl(int64_t* l, double* x, double* jl);
}


void BesselRoots(int64_t nmax, int64_t lmax, double* qln);

double sphericalBesselJ(int64_t l, double x);

void write_qln(char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln);

int load_qln(const char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln);


#endif
