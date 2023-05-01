/*
 * Point.cpp
 *
 *  Created on: Mar 30, 2023
 *      Author: Joye Chen
 */

#include "Point.h"

Double PI = MakeDouble(3.141592653, UU); 

Point MakePoint(Double x, Double y) {
    Point temp; 
    temp.x = x; 
    temp.y = y; 
    return temp;
}

Point MakePoint(double x, double y) {
    Point temp; 
    temp.x = MakeDouble(x, UU); 
    temp.y = MakeDouble(y, UU); 
    return temp; 
}

Point MakePoint(Double x, double y) {
    Point temp; 
    temp.x = x;
    temp.y = MakeDouble(y, UU); 
    return temp; 
}

Point MakePoint(double x, Double y) {
    Point temp; 
    temp.x = MakeDouble(x, UU); 
    temp.y = y; 
    return temp; 
}

Point operator +(Point p, Point q) {
    Point temp; 
    temp.x = p.x + q.x; 
    temp.y = p.y + q.y; 
    return temp; 
}

Point operator -(Point p, Point q) {
    Point temp; 
    temp.x = p.x - q.x; 
    temp.y = p.y - q.y;
    return temp; 
}

Point operator *(Double c, Point p) {
    Point temp; 
    temp.x = c*p.x; 
    temp.y = c*p.y; 
    return temp; 
}

Double distSquared(Point p, Point q) {
    return Pow(p.x - q.x, 2) + Pow(p.y - q.y, 2); 
}

// computes one of the intersection points between two circles of radius r centered at p1 and p2
Point intersectTwoCirc1( Point p1, Point p2, double r ) {
    Double distSq = distSquared( p1, p2 );
    if (distSq > pow(2*r, 2)) { 
        printf("Empty intersection. \n");
		exit(1);
    } 
    Point p3; 
    p3.x = 0.5 * (p1.x + p2.x); 
    p3.y = 0.5 * (p1.y + p2.y); 

    Double b = distSquared(p1, p3); 
    double h = sqrt(r*r - b.value);

    Point p4; 
    p4.x = p3.x - h * (p2.y - p1.y) * (1 / (2 * sqrt(b)));
    p4.y = p3.y - h * (p2.x - p1.x) * (1 / (2 * sqrt(b)));
    return p4; 
}

Point intersectTwoCirc2( Point p1, Point p2, double r) {
    Double distSq = distSquared( p1, p2 );
    if (distSq > pow(2*r, 2)) { 
        printf("Empty intersection. \n");
		exit(1);
    }
        Point p3; 
        p3.x = 0.5 * (p1.x + p2.x); 
        p3.y = 0.5 * (p1.y + p2.y); 

        Double b = distSquared(p1, p3); 
        double h = sqrt(r * r - b.value);

        Point p4; 
        p4.x = p3.x + h * (p2.y - p1.y) * (1 / (2 * sqrt(b)));
        p4.y = p3.y + h * (p2.x - p1.x) * (1 / (2 * sqrt(b)));
        return p4; 
}
    
Point intersectTwoCircGen1( Point p1, Point p2, double r1, double r2) {
    Double distSq = distSquared( p1, p2 );
    Double dist = sqrt(distSq); 

    Double x = (r1*r1 - r2*r2 + distSq) / (2 * dist); 
    if (p1.x > p2.x) { x = -x; }
    double alpha = acos(x.value / r1); 
    double beta = atan((p2.y.value - p1.y.value) / (p2.x.value - p1.x.value)); 

    Point temp; 
    temp.x = r1 * cos(alpha + beta) + p1.x;
    temp.y = r1 * sin(alpha + beta) + p1.y; 

    return temp; 
}

Point intersectTwoCircGen2( Point p1, Point p2, double r1, double r2) {
    Double distSq = distSquared( p1, p2 );
    Double dist = sqrt(distSq);

    Double x = (r1*r1 - r2*r2 + distSq) / (2 * dist); 
    if (p1.x > p2.x) { x = -x; } 
    double alpha = acos(x.value / r1); 
    double beta = atan((p2.y.value - p1.y.value) / (p2.x.value - p1.x.value)); 

    // printf("            alpha = %4.2f, beta = %4.2f.\n", alpha, beta); 


    Point temp; 
    temp.x = r1 * cos(- alpha + beta) + p1.x;
    temp.y = r1 * sin(- alpha + beta) + p1.y;
    // printf("            temp = (%4.2f, %4.2f).\n", temp.x.value, temp.y.value); 

    return temp; 
}

// computes intersection point between circle of radius r centered at p and line passing through origin and p
Point intersectCircLine1( Point p, double r) {
    Point temp; 
    Double A = 1 + Pow(p.y, 2) / Pow(p.x, 2);
    Double B = -2 * p.x - 2 * Pow(p.y, 2) / p.x; 
    Double C = Pow(p.x, 2) + Pow(p.y, 2) - r * r;

    temp.x = (-B - sqrt(B*B - 4 * A * C))/(2 * A); 
    temp.y = temp.x * p.y / p.x; 
    return temp; 
}

Point intersectCircLine2( Point p, double r ) {
    Point temp; 
    Double A = 1 + Pow(p.y, 2) / Pow(p.x, 2);
    Double B = -2 * p.x - 2 * Pow(p.y, 2) / p.x; 
    Double C = Pow(p.x, 2) + Pow(p.y, 2) - r * r;

    temp.x = (-B + sqrt(B*B - 4 * A * C))/(2 * A); 
    temp.y = temp.x * p.y / p.x; 
    return temp; 
}

// computes intersection points between circle of radius r centered at p and line y = mx + b
Point intersectCircLinGen1( Point p, double r, Double m, Double b) {
    Point temp; 
    Double A = m * m + 1.0; 
    Double B = -2 * p.x + 2*m*(b - p.y); 
    Double C = p.x * p.x + (b - p.y) * (b - p.y) - r*r;

    temp.x = (-B - sqrt(B.value*B.value - 4 * A.value * C.value))/(2 * A); 
    temp.y = m * temp.x + b; 
    return temp; 
}

Point intersectCircLinGen2( Point p, double r, Double m, Double b) {
    Point temp; 
    Double A = m * m + 1.0;  
    Double B = -2 * p.x + 2*m*(b - p.y); 
    Double C = p.x * p.x + (b - p.y) * (b - p.y) - r*r; 

    temp.x = (-B + sqrt(B.value*B.value - 4 * A.value * C.value))/(2 * A); 
    temp.y = m * temp.x + b; 
    return temp; 
}

bool isIntersectTwoCirc(Point p1, Point p2, double r1, double r2) {
    double distSq = distSquared(p1, p2).value; 
    if (distSq > (r1+r2)*(r1+r2)) {
        return false;
    }
    else { return true; }
}

bool isIntersectCircLin( Point p, double r, Double m, Double b) {
    Double A = m * m + 1.0;  
    Double B = -2 * p.x + 2*m*(b - p.y); 
    Double C = p.x * p.x + (b - p.y) * (b - p.y) - r*r;

    double disc = B.value*B.value - 4*A.value*C.value; 
    // printf("disc = %6.5f.\n", disc); 
    if (disc < 0.0) {
        return false; 
    }

    else { return true; }
}

bool isInP( Double l, Double b, Double h, Point p ) {
    if (p.x < 0.0 + UU && p.y < -h / b * p.x) { return false; }
    if (p.x > 0.0 && p.x < l - b - UU && (p.y > h || p.y < 0)) { return false; }
	if (p.x > l - b - UU && p.y > -h / b * p.x + l * h / b) { return false; }
    else { return true; }
}

Point translateInP( Double l, Double b, Double h, Point p ) {
    if (isInP(l, b, h, p) == true) { return p; } 
    else {
        double temp = p.y.value / h.value; 
        temp = floor(temp); 
        p.y = p.y - temp * h; // p.y should now be between 0.0 and h as desired. 
        p.x = p.x + temp * b; 
        
        double minX = (p.y * (-b / h)).value; 
        p.x = p.x + minX; 
        temp = p.x.value / l.value; 
        temp = floor(temp); 
        p.x = p.x - temp * l - minX; // p.x should now be between minX and minX + l as desired.

        return p; 
    }
}

// returns angle of p relative to o with value in [0, 2PI)
double angle( Point o, Point p ) {
    double dist = sqrt(distSquared(o, p).value); 
    double xrel = (p.x.value - o.x.value)/dist; 
    double yrel = (p.y.value - o.y.value)/dist; 

    if (yrel < 0) {
        // p in third or fourth quadrant
        return 2*PI.value - acos(xrel); 
    }

    else { 
        // p in first or second quadrant
        return acos(xrel);
    }
}