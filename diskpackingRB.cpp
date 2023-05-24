//============================================================================
// Name        : diskpackingnew.cpp
// Author      : Joye Chen
// Version     :
// Copygight   : Your copygight notice
// Descgiption : This program is a rigorized version of findDiam3.cpp using 
//               standard interval arithmetic (see Double.cpp). It shows that
//               o_1^{RB} > D(v, vR, o1^{RG}); combined with diskpackingnew.cpp
//               we can prove Lemma 5.5. 
//============================================================================


#include <iostream>

// #include "Double.cpp"
#include "Double.h"
#include "Point.cpp"
#include "Point.h"


using namespace std;

Double expApprox( Double x ); // polynomial approximation of exp(x) for x near 0.0 
Double logApprox( Double x ); // polynomial approximation of log(x) for x near ???
Point getFund( Point p, Point o, Point m, Point n); // get fundamental parallelogram for \langle m, n \rangle such that 
// vertices are on translates of o && p contained in interior
// assumes p and o are both contained in some fund. par. (but with vertices not at translates of o)
// also assumes m and n are of form m = (l, 0), n = (h, -b)
// useful for checking p (center of a G) is far enough from translates of o (center of an FSG)
void checkArc(Point FSR, Double angle1, Double angle2, Double distR, Double distG, 
    double width, Double l, Double b, Double h, Point G1, Point G2, Point G4);


Double expApprox(Double x) { // 
    int n = 10; // use first ten terms of Taylor expansion about 0.0
    int fact = n; // compute (n+1)!
    Double temp = MakeDouble(1.0, UU); 
    for (int i = n; i > 0; --i) {
        temp = 1 + x * temp / i; 
        fact = fact * i; 
    }
    fact = (n + 1) * fact; 
    // update error, assuming |x| < 1
    temp.error = 1.0 / fact;
    return temp; 
}

Point getFund(Point p, Point o, Point m, Point n) {
    if (p.y < o.y) {
        while (p.y < o.y) { o = o - n; } 
    }
    if (p.x < o.x) {
        while (p.x < o.x) { o = o - m; }
    }
    return o; 
}

// given FSR and angle1, angle2, check all placements of R along an arc distance distR away 
// from FSR between angle1 and angle2
// check that all such R are "close" to translates of the Gi or translates of FSR
void checkIso(Point FSR, Double angle1, Double angle2, Double distR, Double distG, 
    double width, Double l, Double b, Double h, Point G1, Point G2, Point G4) {
    Double angle; 
    Point R;
    for (angle = angle1; angle > angle2 - width; angle = angle - 2.0 * width) {
        R.x = FSR.x + distR * cosApprox(angle);
        R.y = FSR.y + distR * sinApprox(angle);

        // if G does not lie in P, translate so it does
        if (!isInP(l, b, h, R)) { R = translateInP(l, b, h, R); } 
        // check distances between R and Gi and between R and translates of FSR
        if (isIsolated(R, G1, G2, G4, distG)) { // if G is far enough from FSGs
            if (isIsolated(R, getFund(R, FSR, G2, G4), G2, G4, distR)) { 
                // if R is also far enough from FSRs, COUNTEREXAMPLE!
                printf("l = %6.5f, b = %6.5f, h = %6.5f, FSR = (%4.2f, %4.2f), R = (%4.2f, %4.2f).\n", l.value, b.value, h.value, FSR.x.value, FSR.y.value, R.x.value, R.y.value);
                printf("disttoG = %6.5f, disttoR = %6.5f.\n",distG.value, distR.value); 
                printf("getFund = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", getFund(R, FSR, G2, G4).x.value, getFund(R, FSR, G2, G4).y.value, getFund(R, FSR, G2, G4).x.error, getFund(R, FSR, G2, G4).y.error);
                isIsolatedPrint(R, getFund(R, FSR, G2, G4), G2, G4, distR);
                exit(1);
            }
        }
    }
}

void checkArc(Point FSR, Point Gi, Point Gj, Double distR, Double distG, double width, Double l, Double b, Double h) {
    Point R; 
    Double ZERO = MakeDouble(0.0, UU); 
    Point G1 = MakePoint(ZERO, ZERO); 
    Point G2 = MakePoint(l, ZERO); 
    Point G4 = MakePoint(-b, h); 

    Point temp1 = intersectTwoCircGen1(FSR, Gi, distR, distG); 
    Point temp2 = intersectTwoCircGen2(FSR, Gj, distR, distG); 
    Double acos1 = (temp1.x - FSR.x) / distR; 
    
    Double phi1 = acosApprox((temp1.x - FSR.x) / distR); 
    Double phi2 = acosApprox((temp2.x - FSR.x) / distR); 
    if (temp1.y - FSR.y < 0) { phi1 = -phi1; }
    if (temp2.y - FSR.y < 0) { phi2 = -phi2; }

    if (phi2 > phi1 && abs(phi1 - phi2) > 0.01 ) {
        R.x = FSR.x + distR * cosApprox(0.5 * (phi1 + phi2)); 
        R.y = FSR.y + distR * sinApprox(0.5 * (phi1 + phi2)); 
        if (distSquared(R, Gi) > distR * distR) {
            phi2 = phi2 - 2 * PI; 
        } 
    }

    // pgintf("    phi1 = %4.2f +- %6.5f, phi2 = %4.2f +- %6.5f.\n", phi1.value, phi1.error, phi2.value, phi2.error);
    if (phi2 < phi1 || abs(phi1 - phi2) > 0.01) { checkIso(FSR, phi1, phi2, distG, distR, width, l, b, h, G1, G2, G4); }
}

int main() {
    // define parameters for each choice of parallelogram P
    Double l, h, b; // resp. length, height, horizontal bias
    Double alpha; // angle \in [0, \pi/2] of P, = atan(h/b)
    Point FSR, R; // centers of full-sized red horoball = FSR and second-order red horoball = R
    Point G1, G2, G3, G4; // centers of FSGs at vertices of P

    double width = 1.0/512.0; // temporary width setting to reduce runtime
	l.error = width;
	h.error = width;
	b.error = width;

    // constants
    Double ZERO = MakeDouble(0.0, UU); 
    Double ONE = MakeDouble(1.0, UU); 
    Double PI = MakeDouble(3.14159265358979323846, UU);
    Double sqrt2 = MakeDouble(1.4142135623730951, UU);
    Double volM = MakeDouble(5.3335, UU); // max vol(M)
    Double Boroczky = MakeDouble(0.853, UU);
    Double GHMTY = MakeDouble(2.62, UU); 

    Double v; // v = vol(B) = vol(G)
    double vRmin, vRmax; // vR \in [vRmin, vRmax]
    Double Dmin, Dmax; // Dmin = D(v, vRmax, 0), Dmax = D(v, vRmin, 0)


    // this array stores values of v subdividing the interval [1.41, 1.52]
    double vRB[12] = {1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.51, 1.52};
    // this array stores values of vR we want to test for showing o_1^{RB} > D:
    double vRintervals[12][6] = {
        {1.41, 1.47, 1.53, 1.59, 1.65, 1.73},
        {1.42, 1.47, 1.53, 1.59, 1.65, 1.71},
        {1.43, 1.49, 1.55, 1.61, 1.69, -1}, 
        {1.44, 1.49, 1.55, 1.61, 1.67, -1}, 
        {1.45, 1.51, 1.57, 1.65, -1, -1}, 
        {1.46, 1.51, 1.57, 1.63, -1, -1},
        {1.47, 1.53, 1.61, -1, -1, -1},
        {1.48, 1.53, 1.59, -1, -1, -1},
        {1.49, 1.57, -1, -1, -1, -1},
        {1.50, 1.55, -1, -1, -1, -1},
        {1.51, 1.53, -1, -1, -1, -1},
        {1.52, -1, -1, -1, -1, -1}
    };
    // this array stores D(v, vR = (5.3335)(0.853)-2v, 0.0); is used for showing o_1^{RB} > D only:
    // where v = vRB[i][j]
    double DvalRB[12][6] = {
        {0.341267, 0.305054, 0.274275, 0.247205, 0.22287, 0.193648}, // v = 1.41
        {0.322323, 0.294298, 0.264898, 0.238825, 0.215253, 0.193648}, // v = 1.42
        {0.305054, 0.274275, 0.247205, 0.22287, 0.193648, -1}, // v = 1.43
        {0.289118, 0.264898, 0.238825, 0.215253, 0.193648, -1}, // v = 1.44
        {0.274275, 0.264898, 0.247205, 0.193648, -1, -1}, // v = 1.45
        {0.260349, 0.238825, 0.215253, 0.193648, -1, -1}, // v = 1.46
        {0.247205, 0.22287, 0.193648, -1, -1, -1}, // v = 1.47
        {0.23474, 0.215253, 0.193648, -1, -1, -1}, // v = 1.48
        {0.22287, 0.193648, -1, -1, -1, -1}, // v = 1.49
        {0.211528, 0.193648, -1, -1, -1, -1},  // v = 1.50
        {0.200657, 0.193648, -1, -1, -1, -1}, // v = 1.51
        {0.19021, -1, -1, -1, -1, -1} // v = 1.52
    };

    for (int i = 3; i < 12; i++) {
        v = MakeDouble(vRB[i] + 0.005, 0.005); 
        printf("v in [%4.2f, %4.2f].\n", smallDsize(v), bigDsize(v));
        for (int j = 0; j < 6; j++) {
            vRmin = vRintervals[i][j];
            Dmax = MakeDouble(DvalRB[i][j]*1.1, UU);
            if (vRintervals[i][j + 1] < 0) { continue; }
            if (j == 5) { continue; }
            
            vRmax = vRintervals[i][j + 1];
            Dmin = MakeDouble(DvalRB[i][j + 1]*1.1, UU);
        

            Double diamR = expApprox(-Dmax); // min. Euc. diam. of second-order red horoball R
            // constraints on distances between horoballs
            Double distRtoFSR = sqrt(GHMTY * diamR/vRmax); // min. Euc. distance between R and any FSR
            Double distFSGtoFSG = sqrt(GHMTY / v); // min. Euc. distance between FSGs
            Double distRtoFSG = sqrt(diamR); // min. Euc. distance from R to FSG

            printf("    vR in [%4.2f, %4.2f].\n", vRmin, vRmax);
            printf("        distRtoFSR = %6.5f +- %6.5f.\n", distRtoFSR.value, distRtoFSR.error);
            printf("        distRtoFSG = %6.5f +- %6.5f.\n", distRtoFSG.value, distRtoFSG.error);

            Double minLength = sqrt(2.0 * v); // can assume length >= height
            if (distFSGtoFSG > minLength) {
                minLength = distFSGtoFSG;
            }
            Double tempC = (GHMTY / 2.0)/vRmin; 
            // maximum length of parallelogram P we have to check given vol(B) = vol(G) = v
            Double maxLength = sqrt(sqrt(4*Pow(v, 4) / Pow(tempC, 2) + 2 * Pow(v, 2) + Pow(tempC, 2)/4));
            Double minB, maxB; 
            Double maxAlpha; 

            // miscellaneous
            double time; // to store clock time
            Double theta, phi; // increment for choosing placement of FSR, R resp. 
            Point temp1, temp2; // store intersection points (for choosing placement of R)
            Point a1, a2; // used in choosing which intersection points to consider
            Double phi1, phi2; // store angles made by temp1, temp2 relative to center of FSR (for choosing placement of R)

            // boolean values tracking whether Circ(FSR, distRR) intersect Gi (for placing R)
            bool intersectG1; 
            bool intersectG2; 
            bool intersectG3; 
            bool intersectG4;

            // first two layers parametrize P
            // l < maxLength
            for (l = minLength + width; l < maxLength; l = l + 2 * width) {
                h = 2 * v / l; // suffices to check parallelogram of area 2v 
                // printf("minLength = %6.5f +- %6.5f.\n", minLength.value, minLength.error);
                // printf("height = %6.5f+- %6.5f.\n", h.value, h.error);
                if (h < 0.95*minLength) {
                    // printf("h < minLength, (minLength + h) * (minLength - h) = %6.5f +- %6.5f.\n", ((minLength + h) * (minLength - h)).value, ((minLength + h) * (minLength - h)).error);
                    minB.value = smallDsize(sqrt((minLength + h) * (minLength - h)));
                    minB.error = width;
                }
                else { minB = ZERO; }
                // printf("minB = %6.5f +- %6.5f.\n", minB.value, minB.error);
                maxB = l - sqrt(bigDsize(l) * bigDsize(l) - Pow(2 * smallDsize(v)/l, 2)) + width; // upper bound on horizontal bias

                // check on runtime:
                time = clock();
                printf("    length: %6.5f (time elapsed: %4.2f seconds)\n", l.value, time / CLOCKS_PER_SEC);

                // b < maxB
                for (b = minB + width; b < maxB; b = b + 2 * width) {
                    // check Euc. dists btwn FSRs are >= distRR
                    Double sideLength = sqrt(Pow(b, 2) + Pow(h, 2));
                    if (sideLength < distFSGtoFSG) { continue; }

                    // set centers of FSRs
                    G1 = MakePoint(ZERO, ZERO); 
                    G2 = MakePoint(l, ZERO); 
                    G3 = MakePoint(l - b, h); 
                    G4 = MakePoint(-b, h); 
                    // printf("R4 = (%6.5f, %6.5f).\n", R4.x.value, R4.y.value);

                    // next layer selects for center of FSG; by symmetry can choose FSG to lie on boundary of R4 or on boundary of R1
                    if (isIntersectTwoCirc(G1, G4, ONE, ONE)) {
                        // let intR1R4 be the point of intersection contained in P
                        Point intR1R4 = intersectTwoCircGen1(G1, G4, ONE, ONE);
                        maxAlpha = acosApprox((b + intR1R4.x));
                    }
                    else { maxAlpha = acosApprox(b / sideLength); } 

                    // first choose FSG to lie on boundary(R4)
                    for (theta = ZERO + width; theta < maxAlpha; theta = theta + 2.0 * width) {
                        FSR.x = G4.x + cosApprox(theta);
                        FSR.y = G4.y - sinApprox(theta);

                        // check that FSR is not too close to the Gi
                        if (distSquared(FSR, G2) < 1.0 || distSquared(FSR, G3) < 1.0 || distSquared(FSR, G1) < 1.0) { continue; }

                        // record whether Circ(FSG, distGtoFSG) intersects Circ(Ri, distGtoFSR) for placement of G
                        intersectG1 = isIntersectTwoCirc(FSR, G1, distRtoFSR, distRtoFSG);
                        intersectG2 = isIntersectTwoCirc(FSR, G2, distRtoFSR, distRtoFSG);
                        intersectG3 = isIntersectTwoCirc(FSR, G3, distRtoFSR, distRtoFSG);

                        if (intersectG3) { checkArc(FSR, G4, G3, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else if (intersectG2) { checkArc(FSR, G4, G2, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else if (intersectG1) { checkArc(FSR, G4, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else { checkArc(FSR, G4, G4, distRtoFSR, distRtoFSG, width, l, b, h); }

                        if (intersectG3) {
                            if (intersectG2) { checkArc(FSR, G3, G2, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else if (intersectG1) { checkArc(FSR, G3, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else { checkArc(FSR, G3, G4, distRtoFSR, distRtoFSG, width, l, b, h); }
                        }

                        if (intersectG2) {
                            if (intersectG1) { checkArc(FSR, G2, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else { checkArc(FSR, G2, G4, distRtoFSR, distRtoFSG, width, l, b, h); }
                        }

                        if (intersectG1) { checkArc(FSR, G1, G4, distRtoFSR, distRtoFSG, width, l, b, h); }

                    }

                    maxAlpha = PI - maxAlpha; 
                    for (theta = ZERO + width; theta < maxAlpha; theta = theta + 2.0 * width) {
                        FSR.x= G1.x + cosApprox(theta); 
                        FSR.y = G1.y + sinApprox(theta);

                        // check that FSR is not too close to the Gi
                        if (distSquared(FSR, G2) < 1.0 || distSquared(FSR, G3) < 1.0 || distSquared(FSR, G4) < 1.0) { continue; }

                        // record whether Circ(FSG, distGtoFSG) intersects Circ(Ri, distGtoFSR) for placement of G
                        intersectG4 = isIntersectTwoCirc(FSR, G4, distRtoFSR, distRtoFSG);
                        intersectG2 = isIntersectTwoCirc(FSR, G2, distRtoFSR, distRtoFSG);
                        intersectG3 = isIntersectTwoCirc(FSR, G3, distRtoFSR, distRtoFSG);

                        if (intersectG4) { checkArc(FSR, G1, G4, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else if (intersectG3) { checkArc(FSR, G1, G3, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else if (intersectG2) { checkArc(FSR, G1, G2, distRtoFSR, distRtoFSG, width, l, b, h); }
                        else { checkArc(FSR, G1, G1, distRtoFSR, distRtoFSG, width, l, b, h); }

                        if (intersectG4) {
                            if (intersectG3) { checkArc(FSR, G4, G3, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else if (intersectG2) { checkArc(FSR, G4, G2, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else { checkArc(FSR, G4, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                        }

                        if (intersectG3) {
                            if (intersectG2) { checkArc(FSR, G3, G2, distRtoFSR, distRtoFSG, width, l, b, h); }
                            else { checkArc(FSR, G3, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                        }
                        
                        if (intersectG2) { checkArc(FSR, G2, G1, distRtoFSR, distRtoFSG, width, l, b, h); }
                    }
                }
            }
            
        }
    }
}