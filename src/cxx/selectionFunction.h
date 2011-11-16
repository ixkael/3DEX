/*
 *  selectionFunction.h
 *  MRS3D
 *
 *  Created by François Lanusse on 17/06/11.
 *  Copyright 2011 CEA. All rights reserved.
 *
 */


#ifndef SELECTION_FUNCTION_H
#define SELECTION_FUNCTION_H

#include <math.h>

template<typename T> class selectionFunction{
	
public:
	
	virtual T evaluate(T rho){
		return 1;
	}
};

template<typename T> class gaussianSelection : public selectionFunction<T>{
	
	
public:
	gaussianSelection(T rho0) : 
	rho0_(rho0)
	{}
	
	T evaluate(T rho){
		return 4.0/(sqrt(M_PI)*rho0_*rho0_*rho0_) * exp(-pow(rho/rho0_, 2));
	}
	
private:
	T rho0_;
};

#endif