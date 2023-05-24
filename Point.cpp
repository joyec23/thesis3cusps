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
Point intersectTwoCirc1( Point p1, Point p2, Double r ) {
    Double distSq = distSquared( p1, p2 );
    if (smallDsize(distSq) > bigDsize(Pow(2*r, 2))) { 
        printf("Empty intersection. \n");
		exit(1);
    } 
    Point p3; 
    p3.x = 0.5 * (p1.x + p2.x); 
    p3.y = 0.5 * (p1.y + p2.y); 

    Double b = distSquared(p1, p3); 
    double h = sqrt(r*r - b);

    Point p4; 
    p4.x = p3.x - h * (p2.y - p1.y) * (1 / (2 * sqrt(b)));
    p4.y = p3.y - h * (p2.x - p1.x) * (1 / (2 * sqrt(b)));
    return p4; 
}

Point intersectTwoCirc2( Point p1, Point p2, Double r) {
    Double distSq = distSquared( p1, p2 );
    if (smallDsize(distSq) > bigDsize(Pow(2*r, 2))) { 
        printf("Empty intersection. \n");
		exit(1);
    }
        Point p3; 
        p3.x = 0.5 * (p1.x + p2.x); 
        p3.y = 0.5 * (p1.y + p2.y); 

        Double b = distSquared(p1, p3); 
        double h = sqrt(r * r - b);

        Point p4; 
        p4.x = p3.x + h * (p2.y - p1.y) * (1 / (2 * sqrt(b)));
        p4.y = p3.y + h * (p2.x - p1.x) * (1 / (2 * sqrt(b)));
        return p4; 
}
    
// WARNING: check for non-empty intersection before calling
Point intersectTwoCircGen1( Point p1, Point p2, Double r1, Double r2) {
    Double d = sqrt(distSquared( p1, p2 )); 
    Double l = (r1 * r1 - r2 * r2 + d * d) / (2*d); 
    Double h = sqrt(r1 * r1 - l * l); 

    Double x = (l / d) * (p2.x - p1.x) + (h / d) * (p2.y - p1.y) + p1.x; 
    Double y = (l / d) * (p2.y - p1.y) - (h / d) * (p2.x - p1.x) + p1.y; 
    return MakePoint(x, y); 
}

// WARNING: check for non-empty intersection before calling
Point intersectTwoCircGen2( Point p1, Point p2, Double r1, Double r2) {
    Double d = sqrt(distSquared( p1, p2 )); 
    Double l = (r1 * r1 - r2 * r2 + d * d) / (2*d); 
    Double h = sqrt(r1 * r1 - l * l); 

    Double x = (l / d) * (p2.x - p1.x) - (h / d) * (p2.y - p1.y) + p1.x; 
    Double y = (l / d) * (p2.y - p1.y) + (h / d) * (p2.x - p1.x) + p1.y; 
    return MakePoint(x, y); 
}

// computes intersection point between circle of radius r centered at p and line passing through origin and p
Point intersectCircLine1( Point p, Double r) {
    Point temp; 
    Double A = 1 + Pow(p.y, 2) / Pow(p.x, 2);
    Double B = -2 * p.x - 2 * Pow(p.y, 2) / p.x; 
    Double C = Pow(p.x, 2) + Pow(p.y, 2) - r * r;

    temp.x = (-B - sqrt(B*B - 4 * A * C))/(2 * A); 
    temp.y = temp.x * p.y / p.x; 
    return temp; 
}

Point intersectCircLine2( Point p, Double r ) {
    Point temp; 
    Double A = 1 + Pow(p.y, 2) / Pow(p.x, 2);
    Double B = -2 * p.x - 2 * Pow(p.y, 2) / p.x; 
    Double C = Pow(p.x, 2) + Pow(p.y, 2) - r * r;

    temp.x = (-B + sqrt(B*B - 4 * A * C))/(2 * A); 
    temp.y = temp.x * p.y / p.x; 
    return temp; 
}

// computes intersection points between circle of radius r centered at p and line y = mx + b
Point intersectCircLinGen1( Point p, Double r, Double m, Double b) {
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

bool isIntersectTwoCirc(Point p1, Point p2, Double r1, Double r2) {
    Double distSq = distSquared(p1, p2); 
    if (bigDsize(distSq) > smallDsize((r1+r2)*(r1+r2))) {
        return false;
    }
    Double d = sqrt(distSquared( p1, p2 )); 
    Double l = (r1 * r1 - r2 * r2 + d * d) / (2*d); 
    Double hSq = r1 * r1 - l * l;
    if (smallDsize(hSq) < 0) { return false; }
    else { return true; }
}

bool isIntersectCircLin( Point p, Double r, Double m, Double b) {
    Double A = m * m + 1.0;  
    Double B = -2 * p.x + 2*m*(b - p.y); 
    Double C = p.x * p.x + (b - p.y) * (b - p.y) - r*r;

    Double disc = B * B - 4 * A * C; 
    // printf("disc = %6.5f.\n", disc); 
    if (smallDsize(disc) < 0.0) {
        return false; 
    }
    else { return true; }
}

bool isInP( Double l, Double b, Double h, Point p ) {
    if (p.x < -b || p.x > l) { return false; }
    if (p.y < 0.0 || p.y > h) { return false; }
    if (p.x < 0.0 + UU && p.y < -h / b * p.x) { return false; }
    if ((p.x > -b && p.x < 0.0 + UU) && p.y > h) { return false; }
    if (p.x > 0.0 && p.x < l - b - UU && (p.y > h || p.y < 0)) { return false; }
	if (p.x > l - b && (p.y > -h / b * (p.x - l) || p.y < 0)) { return false; }
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
        p.x = p.x - minX; 
        temp = p.x.value / l.value; 
        temp = floor(temp); 
        p.x = p.x - temp * l + minX; // p.x should now be between minX and minX + l as desired.

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

bool isIsolated(Point p, Point o, Point m, Point n, Double dist) {
    // printf("boop20.\n"); 
    Double distSq = dist * dist; 
    if (distSquared(p, o) < distSq) { // printf("boop21.\n");  
        return false; 
    }
    if (distSquared(p, m + o) < distSq) { // printf("boop22.\n");
        return false; 
    } 
    if (distSquared(p, n + o) < distSq) { // printf("boop23.\n"); 
        return false; 
    }
    if (distSquared(p, m + n + o) < distSq) { 
        // Double distTemp = distSquared(p, m + n + o); 
        // printf("p = (%6.5f, %6.5f) +- (%6.5f, %6.5f), m + n + o = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", p.x.value, p.y.value, p.x.error, p.y.error, (m + n + o).x.value, (m + n + o).y.value, (m + n + o).x.error, (m + n + o).y.error);
        // printf("    distSq = %6.5f +- %6.5f.\n", distTemp.value, distTemp.error); 
        // printf("boop24.\n");
        return false; 
    }
    if (distSquared(p, m + m + n + o) < distSq) {
        return false;
    }
    else { return true; }
}

void isIsolatedPrint(Point p, Point o, Point m, Point n, Double dist) {
    Double distSq = dist * dist; 
    printf("G = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", p.x.value, p.y.value, p.x.error, p.y.error);
    printf("o = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", o.x.value, o.y.value, o.x.error, o.y.error);
    printf("m = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", m.x.value, m.y.value, m.x.error, m.y.error);
    printf("n = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", n.x.value, n.y.value, n.x.error, n.y.error);
    printf("    distSquared(p, o) = %6.5f +- %6.5f. \n", distSquared(p, o).value, distSquared(p, o).error); 
    printf("    distSquared(p, m + o) = %6.5f +- %6.5f.\n", distSquared(p, m + o).value, distSquared(p, m + o).error); 
    printf("    distSquared(p, n + o) = %6.5f +- %6.5f.\n", distSquared(p, n + o).value, distSquared(p, n + o).error); 
    printf("    distSquared(p, m + n + o) = %6.5f +- %6.5f.\n", distSquared(p, m + n + o).value, distSquared(p, m + n + o).error);

}
