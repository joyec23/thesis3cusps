/*
 * Double.h
 *
 *  Created on: Mar 11, 2023
 *      Author: Marc Lackenby & Rob Meyerhoff
 */

#ifndef __myDouble_h__
#define __myDouble_h__
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#define UU  (DBL_EPSILON/2)
typedef struct {
	    double value;
	    double error;
	}   Double;

extern	Double 		MakeDouble		(double x, double y);
extern	Double 		operator +		(Double x, Double y);
extern	Double 		operator +		(double x, Double y);
extern	Double 		operator +		(Double x, double y);
extern	Double	 	operator -		(Double x, Double y);
extern	Double	 	operator -		(double x, Double y);
extern	Double	 	operator -		(Double x, double y);
extern	Double	 	operator -		(Double x);
extern	Double		operator *		(Double x, Double y);
extern	Double		operator *		(double x, Double y);
extern	Double		operator *		(Double x, double y);
extern	Double	 	operator /		(Double x, Double y);
extern	Double 		operator /		(Double x, double y);
extern	Double 		operator /		(double x, Double y);
extern	Double		sqrt			(Double x);
extern	Double 		abs				(Double x);
extern	Double 		hypot			(Double x, Double y);
extern	double 		hypot			(double x, double y);
extern  Double 		expo			(double x, int j);
extern  Double 		Pow				(Double x, int j);
extern	double 		bigDsize		(Double x);
extern	double 		smallDsize		(Double x);
extern	int 		operator <		(Double x, Double y);
extern	int			operator <		(Double x, int y);
extern	int			operator <		(Double x, double y);
extern	int 		operator >		(Double x, Double y);
extern	int			operator >		(Double x, int y);
extern	int			operator >		(Double x, double y);
extern	void 		printDouble		(Double x);

extern	Double		shift			(Double x); 
extern	Double		cosApprox		(Double x);
extern	Double		sinApprox		(Double x); 
extern	Double		acosApprox		(Double x);
extern 	Double		asinApprox		(Double x);
extern 	Double		atanApprox		(Double x); 


#endif
