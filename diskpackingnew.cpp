//============================================================================
// Name        : diskpackingnew.cpp
// Author      : Joye Chen
// Version     :
// Copyright   : Your copyright notice
// Description : This program is a rigorized version of findDiam3.cpp using 
//               standard interval arithmetic (see Double.cpp). It shows that
//               o_1^{BG} > D(v, vR, o1^{RG}); see the proof of Lemma 5.5.
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
void checkArc(Point FSG, Double angle1, Double angle2, Double distG, Double distR, 
    double width, Double l, Double b, Double h, Point R1, Point R2, Point R4);


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

// given FSG and angle1, angle2, check all placements of G along an arc distance distG away 
// from FSG between angle1 and angle2
// check that all such G are "close" to translates of the Ri or translates of FSG
void checkIso(Point FSG, Double angle1, Double angle2, Double distG, Double distR, 
    double width, Double l, Double b, Double h, Point R1, Point R2, Point R4) {
    Double angle; 
    Point G;
    for (angle = angle1; angle > angle2 - width; angle = angle - 2.0 * width) {
        G.x = FSG.x + distG * cosApprox(angle);
        G.y = FSG.y + distG * sinApprox(angle);

        // if G does not lie in P, translate so it does
        if (!isInP(l, b, h, G)) { G = translateInP(l, b, h, G); } 
        // check distances between G and Ri and between G and translates of FSG
        if (isIsolated(G, R1, R2, R4, distR)) { // if G is far enough from FSRs
            if (isIsolated(G, getFund(G, FSG, R2, R4), R2, R4, distG)) { 
                // if G is also far enough from FSGs, COUNTEREXAMPLE!
                printf("l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                printf("disttoR = %6.5f, disttoG = %6.5f.\n",distR.value, distG.value); 
                printf("getFund = (%6.5f, %6.5f) +- (%6.5f, %6.5f).\n", getFund(G, FSG, R2, R4).x.value, getFund(G, FSG, R2, R4).y.value, getFund(G, FSG, R2, R4).x.error, getFund(G, FSG, R2, R4).y.error);
                isIsolatedPrint(G, getFund(G, FSG, R2, R4), R2, R4, distG);
                exit(1);
            }
        }
    }
}

void checkArc(Point FSG, Point Ri, Point Rj, Double distG, Double distR, double width, Double l, Double b, Double h) {
    Point G; 
    Double ZERO = MakeDouble(0.0, UU); 
    Point R1 = MakePoint(ZERO, ZERO); 
    Point R2 = MakePoint(l, ZERO); 
    Point R4 = MakePoint(-b, h); 

    Point temp1 = intersectTwoCircGen1(FSG, Ri, distG, distR); 
    // printf("checkpoint000.\n");
    Point temp2 = intersectTwoCircGen2(FSG, Rj, distG, distR); 
    // printf("checkpoint0.\n");
    Double acos1 = (temp1.x - FSG.x) / distG; 
    // printf("acos1 input = %6.5f +- %6.5f.\n", acos1.value, acos1.error);
    
    Double phi1 = acosApprox((temp1.x - FSG.x) / distG); 
    Double phi2 = acosApprox((temp2.x - FSG.x) / distG); 
    if (temp1.y - FSG.y < 0) { phi1 = -phi1; }
    if (temp2.y - FSG.y < 0) { phi2 = -phi2; }
    // printf("checkpoint1.\n");

    if (phi2 > phi1 && abs(phi1 - phi2) > 0.01 ) {
        // printf("checkpoint2.\n");
        G.x = FSG.x + distG * cosApprox(0.5 * (phi1 + phi2)); 
        G.y = FSG.y + distG * sinApprox(0.5 * (phi1 + phi2)); 
        if (distSquared(G, Ri) > distG * distG) {
            phi2 = phi2 - 2 * PI; 
        } 
    }

    // printf("    phi1 = %4.2f +- %6.5f, phi2 = %4.2f +- %6.5f.\n", phi1.value, phi1.error, phi2.value, phi2.error);
    if (phi2 < phi1 || abs(phi1 - phi2) > 0.01) { checkIso(FSG, phi1, phi2, distG, distR, width, l, b, h, R1, R2, R4); }
}


int main() {
    // define parameters for each choice of parallelogram P
    Double l, h, b; // resp. length, height, horizontal bias
    Double alpha; // angle \in [0, \pi/2] of P, = atan(h/b)
    Point FSG, G; // centers of full-sized green horoball = FSG and second-order green horoball = G
    Point R1, R2, R3, R4; // centers of FSRs at vertices of P

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

    // idea is to eliminate possibility of valid configs for v \in [v[i], v[i + 1]]

    Double v; // v = vol(B) = vol(G)
    Double vR; // v_R = vol(R)
    Double D; // D(v, v_R, o_1^{RG})
    double volWidth = 0.005; 

    // this array stores values of v subdividing the interval [1.41, 1.52]
    double vBG[12] = {1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.51, 1.52};
    // this array stores values of D(v, vR = v, 0.0):
    double Dval[12] = {0.3413, 0.3223, 0.3051, 0.2891, 0.2743, 0.2603, 0.2472, 0.2347, 0.2229, 0.2115, 0.2007, 0.1902};


    for (int i = 0; i < 12; i++ ) {
        v = MakeDouble(vBG[i] + 0.005, 0.005); 
        vR = volM * Boroczky - 2*v; 
        D = MakeDouble(Dval[i] * 1.2, 0.001); // error term added to D to account for potential floating-point error
        Double diamG = expApprox(-D); // Euclidean height of G
        // constraints on distances between horoballs
        Double distGtoFSGsq = GHMTY * diamG / v; // min Euc. distance squared from G to FSG 
        Double distGtoFSG = sqrt(distGtoFSGsq); 
        Double distFSRtoFSR= sqrt(GHMTY / vR); // min Euc. distance between FSRs
        Double distGtoFSR = sqrt(diamG); // min Euc. distance from G to FSR 

        printf("distGtoFSG = %6.5f +- %6.5f.\n", distGtoFSG.value, distGtoFSG.error);
        printf("distGtoFSR = %6.5f +- %6.5f.\n", distGtoFSR.value, distGtoFSR.error);

        Double minLength = sqrt(2.0 * v); // can assume length >= height
            if (distFSRtoFSR > minLength) {
                minLength = distFSRtoFSR;
            }
        Double tempC = (GHMTY / 2.0)/vR; 
        // maximum length of parallelogram P we have to check given vol(B) = vol(G) = v
        Double maxLength = sqrt(sqrt(4*Pow(v, 4) / Pow(tempC, 2) + 2 * Pow(v, 2) + Pow(tempC, 2)/4));
        Double minB, maxB; 
        Double maxAlpha;

        // miscellaneous
        double time; // to store clock time
        Double theta, phi; // increment for choosing placement of FSG, G resp. 
        Point temp1, temp2; // store intersection points (for choosing placement of G)
        Point a1, a2; // used in choosing which intersection points to consider
        Double phi1, phi2; // store angles made by temp1, temp2 relative to center of FSG (for choosing placement of G)
        int counterFSG, counterG; // to sample FSGs and Gs for printing


        // boolean values tracking whether Circ(FSG, distGG) intersect Ri (for placing G)
        bool intersectR1; 
        bool intersectR2; 
        bool intersectR3; 
        bool intersectR4;

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
            printf("length: %6.5f (time elapsed: %4.2f seconds)\n", l.value, time / CLOCKS_PER_SEC);

            // b < maxB
            for (b = minB + width; b < maxB; b = b + 2 * width) {
                // check Euc. dists btwn FSRs are >= distRR
                Double sideLength = sqrt(Pow(b, 2) + Pow(h, 2));
                if (sideLength < distFSRtoFSR) { continue; }

                counterFSG = 0; 

                // set centers of FSRs
                R1 = MakePoint(ZERO, ZERO); 
                R2 = MakePoint(l, ZERO); 
                R3 = MakePoint(l - b, h); 
                R4 = MakePoint(-b, h); 
                // printf("R4 = (%6.5f, %6.5f).\n", R4.x.value, R4.y.value);

                // next layer selects for center of FSG; by symmetry can choose FSG to lie on boundary of R4 or on boundary of R1
                if (isIntersectTwoCirc(R1, R4, ONE, ONE)) {
                    // let intR1R4 be the point of intersection contained in P
                    Point intR1R4 = intersectTwoCircGen1(R1, R4, ONE, ONE);
                    maxAlpha = acosApprox((b + intR1R4.x));
                }
                else { maxAlpha = acosApprox(b / sideLength); } 

                // first choose FSG to lie on boundary(R4)
                for (theta = ZERO + width; theta < maxAlpha; theta = theta + 2.0 * width) {
                    counterFSG++;
                    counterG = 0; 
                    FSG.x = R4.x + cosApprox(theta);
                    FSG.y = R4.y - sinApprox(theta);
                    // if (counterFSG % 10 == 0) { printf("    FSG = (%6.5f, %6.5f). \n", FSG.x.value, FSG.y.value); }

                    // check that FSG is not too close to R1, R2, or R3
                    if (distSquared(FSG, R2) < 1.0 || distSquared(FSG, R3) < 1.0 || distSquared(FSG, R1) < 1.0) { continue; }

                    // record whether Circ(FSG, distGtoFSG) intersects Circ(Ri, distGtoFSR) for placement of G
                    intersectR1 = isIntersectTwoCirc(FSG, R1, distGtoFSG, distGtoFSR);
                    intersectR2 = isIntersectTwoCirc(FSG, R2, distGtoFSG, distGtoFSR);  
                    intersectR3 = isIntersectTwoCirc(FSG, R3, distGtoFSG, distGtoFSR); 
                    // note intersectR4 = true; 

                    // compute arcs of Circ(FSG, distGG) which lie outside of disks centered at the Ri 
                    // with radius distGtoFSR, working clockwise from R4
                    // for each arc:
                        // determine whether it is covered by the involved Ri, R(i + 1)
                        // determine whether it is covered by translates of Circ(FSG, distGG), working clockwise

                    // start with arcs having an endpoint on R4
                    if (intersectR3) {
                        // printf("checking arc(R4, R3).\n");
                        checkArc(FSG, R4, R3, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    // check G along arc from R4 to R2
                    else if (intersectR2) {
                        // printf("checking arc(R4, R2).\n");
                        checkArc(FSG, R4, R2, distGtoFSG, distGtoFSR, width, l, b, h);
                    }

                    // check G along arc from R4 to R1
                    else if (intersectR1) {
                        // printf("checking arc(R4, R1).\n");
                        checkArc(FSG, R4, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    // check G along arc from R4 to R4 
                    else {
                        // printf("checking arc(R4, R4).\n");
                        checkArc(FSG, R4, R4, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    // next check arcs w/ endpoint on boundary(R3)
                    if (intersectR3) {             
                        // check G along arc from R3 to R2
                        if (intersectR2) {
                            // printf("checking arc(R3, R2).\n");
                            checkArc(FSG, R3, R2, distGtoFSG, distGtoFSR, width, l, b, h); 

                        }

                        // check G along arc from R3 to R1
                        else if (intersectR1) {
                            // printf("checking arc(R3, R1).\n");
                        checkArc(FSG, R3, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }

                        // check G along arc from R3 to R4
                        else {
                            // printf("checking arc(R3, R4).\n");
                            checkArc(FSG, R3, R4, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }
                    }

                    // next check arcs w/ endpoint on boundary(R2)
                    if (intersectR2) {
                        // check G along arc from R2 to R1
                        if (intersectR1) {
                            // printf("checking arc(R2, R1).\n");
                            checkArc(FSG, R2, R1, distGtoFSG, distGtoFSR, width, l, b, h);
                        }

                        // check G along arc from R2 to R4
                        else {
                            // printf("checking arc(R2, R4).\n");
                            checkArc(FSG, R2, R4, distGtoFSG, distGtoFSR, width, l, b, h);
                        }
                    }
                    
                    // finally, check arc w/ endpoints on boundary(R1) and boundary(R4)
                    if (intersectR1) {
                        // printf("checking arc(R2, R1).\n");
                        checkArc(FSG, R1, R4, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                }


                // next choose FSG to lie on boundary(R1)
                maxAlpha = PI - maxAlpha; 
                for (theta = ZERO + width; theta < maxAlpha; theta = theta + 2.0 * width) {
                    counterFSG++; 
                    counterG = 0; 
                    FSG.x= R1.x + cosApprox(theta); 
                    FSG.y = R1.y + sinApprox(theta);
                    // if (counterFSG % 100 == 0) { printf("    FSG = (%6.5f, %6.5f). \n", FSG.x.value, FSG.y.value); }

                    // check that FSG is not too close to R2, R3, or R4
                    if (distSquared(FSG, R2) < 1.0 || distSquared(FSG, R3) < 1.0 || distSquared(FSG, R4) < 1.0) { continue; }

                    // record whether Circ(FSG, distGtoFSG) intersects Circ(Ri, distGtoFSR) for placement of G
                    // note intersectR1 = true;
                    intersectR2 = isIntersectTwoCirc(FSG, R2, distGtoFSG, distGtoFSR);  
                    intersectR3 = isIntersectTwoCirc(FSG, R3, distGtoFSG, distGtoFSR); 
                    intersectR4 = isIntersectTwoCirc(FSG, R4, distGtoFSG, distGtoFSR); 

                    // compute arcs of Circ(FSG, distGG) which lie outside of disks centered at the Ri 
                    // with radius distGtoFSR, working clockwise from R4
                    // for each arc:
                        // determine whether it is covered by the involved Ri, R(i + 1)
                        // determine whether it is covered by translates of Circ(FSG, distGG), working clockwise
                    // start with arcs having an endpoint on R1

                    if (intersectR4) {
                        checkArc(FSG, R1, R4, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }
                    
                    else if (intersectR3) {
                        checkArc(FSG, R1, R3, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    else if (intersectR2) {
                        checkArc(FSG, R1, R2, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    else {
                        checkArc(FSG, R1, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }

                    // check arcs having an endpoint on R4
                    if (intersectR4) {
                        if (intersectR3) {
                            checkArc(FSG, R4, R3, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }

                        else if (intersectR2) {
                            checkArc(FSG, R4, R2, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }

                        else {
                            checkArc(FSG, R4, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }
                    }

                    // check arcs having an endpoint on R3
                    if (intersectR3) {
                        if (intersectR2) {
                            checkArc(FSG, R3, R2, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }

                        else {
                            checkArc(FSG, R3, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                        }
                    }

                    // finally, check arc from R2 to R1
                    if (intersectR2) {
                        checkArc(FSG, R2, R1, distGtoFSG, distGtoFSR, width, l, b, h); 
                    }
                    
                } 
            }  
        }
    }
}