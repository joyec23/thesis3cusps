/*
 * Point.h
 *
 *  Created on: Mar 30, 2023
 *      Author: Joye Chen
 */

#ifndef __myPoint_h__
#define __myPoint_h__
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "Double.cpp"
#include "Double.h"
typedef struct {
	    Double x;
	    Double y;
	}   Point;

extern	Point 		MakePoint		(Double x, Double y);
extern	Point 		MakePoint		(Double x, double y);
extern	Point 		MakePoint		(double x, Double y);
extern	Point 		MakePoint		(double x, double y);
extern	Point 		operator +		(Point p, Point q);
extern  Point       operator -      (Point p, Point q); 
extern  Point       operator *      (Double c, Point p); 
extern  Double      distSquared     (Point p, Point q);

extern 	Point 	intersectTwoCirc1		( Point p1, Point p2, Double r); 
extern	Point 	intersectTwoCirc2		( Point p1, Point p2, Double r); 
extern	Point 	intersectTwoCircGen1	( Point p1, Point p2, Double r1, Double r2); 
extern	Point 	intersectTwoCircGen1	( Point p1, Point p2, Double r1, Double r2); 
extern	Point 	intersectCircLine1		( Point p, Double r ); 
extern	Point 	intersectCircLine2		( Point p, Double r );
extern	Point 	intersectCircLinGen1	( Point p, Double r, Double m, Double b); 
extern	Point 	intersectCircLinGen2	( Point p, Double r, Double m, Double b); 
extern 	bool 	isIntersectTwoCirc		( Point p1, Point p2, double r1, double r2); 
extern	bool 	isIntersectCircLin		( Point p, double r, Double m, Double b); 

extern	bool 	isInP					( Double l, Double b, Double h, Point p); 
extern  Point	translateInP			( Double l, Double b, Double h, Point p); 

extern  double 	angle					( Point o, Point p ); 
extern bool 	isIsolated				( Point p, Point o, Point m, Point n, Double dist);
extern void 	isIsolatedPrint			( Point p, Point o, Point m, Point n, Double dist); // for debugging use


#endif
