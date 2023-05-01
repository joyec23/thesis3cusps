//============================================================================
// Name        : findDiam2.cpp
// Author      : Joye Chen
// Version     :
// Copyright   : Your copyright notice
// Description : This program tests that, subject to the constraint
// 					vol(M) <= 5.3335 and vol(blue) = cuspVol, there are no valid
//					configurations of an FSR, FSG, and second-order horoball G
//					with diameter e^{-D}. It improves upon findDiam1.cpp by testing 
//                  positions of G lying along the circle of radius distGG centered
//                  at FSG. 
//============================================================================


#include <iostream>
// #include "Double.cpp"
#include "Double.h"
#include "Point.cpp"
#include "Point.h"

using namespace std;

int main() { // the first part of the program is the same as diskpacking.cpp (up through defining P)
	// l = length, h = height, b = horiz skew
	// (x_FSG, y_FSG) = coords of full-sized green, (x_MG, y_MG) = coords of microgreen
	double ll, bb; 
	Double l, h, b, minSL, minLength, SLsquared;
    Point FSG, G; 
    Point R1, R2, R3, R4; 

	double width = 1.0/512.0; // temporary width setting to reduce runtime
	l.error = width;
	h.error = width;
	b.error = width;

    // constants
    Double ZERO = MakeDouble(0.0, UU); 
    Double PI; PI.value = 3.14159265358979323846; PI.error = UU;
    // Double Boroczky; Boroczky.value = 0.853; Boroczky.error = UU;
    // Double sqrt2; sqrt2.value = 1.4142135623730951; sqrt2.error = UU;

	Double cuspVol = MakeDouble(1.41, UU); // possible blue, green cusp volume
    Double cuspVolR = MakeDouble(1.721, UU); // possible red cusp volume
	Double diamMG = MakeDouble(0.8, UU);
	Double distGGsq = diamMG * 2.62 / cuspVolR; // minimum possible distsquared from MG to FSG (or MR to FSR)
    Double distGG = sqrt(distGGsq); 
    printf("distGG = %6.5f.\n", sqrt(distGGsq.value)); 
	
    Double distRGsq = diamMG; // minimum possible dist squared from MG to the FSRs
    Double distRG = sqrt(diamMG); // min possible dist from MG to the FSRs
    printf("distRG = %6.5f.\n", distRG.value); 
    
    double alpha, alpha1, alpha2; // alpha is angle of parallelogram, angle1, angle2 for use in selecting position of G later

    Double diameter = ZERO; // diameter of valid region (for placing FSG and G)
    Double diamTemp; 
	bool existsConfig;
	double time;
    
    Double theta; // later used as an increment to choose FSG 
    Double phi; // later used as an increment to choose G
    theta.error = width;
    phi.error = width; 

    // boolean values tracking whether Circ(FSG, distGG) intersect Ri
    bool intersectR1; 
    bool intersectR2; 
    bool intersectR3; 
    bool intersectR4;  
    Point temp; Point temp1; Point temp2; // for later
    Double tempDist; // for later

	minLength = sqrt(2.0 * cuspVol); // min length, approx 1.7415
	minSL = sqrt(2.62 / cuspVol); // min side length, approx 1.3144 for cuspVol = 1.5165

    int counter; 

	// first two layers parametrize P
    // should start from sqrt(2*cuspVol.value) + width + 0.0/256.0
	for (ll = sqrt(2*cuspVol.value) + width; ll < smallDsize(minLength) + 160.0/256.0; ll = ll + 2.0 * width) {
		l.value = ll;
		time = clock();
		printf("length: %6.5f (time elapsed: %4.2f seconds)\n", l.value, time / CLOCKS_PER_SEC);

		h = 2 * cuspVol / l; 

		for (bb = 0.1 + width; bb < min(0.5 * ll, ll - sqrt(ll * ll - h.value * h.value)) + width; bb = bb + 2.0 * width) {
			time = clock();
			// printf("    horiz bias: %6.5f (time elapsed: %4.2f seconds)\n", bb, time / CLOCKS_PER_SEC);
			b.value = bb;
            time = clock();
            diameter = ZERO; // reset diameter for each new parallelogram
            existsConfig = false; 
            // printf("    b =  %6.5f (time elapsed: %4.2f seconds)\n", b.value, time / CLOCKS_PER_SEC);

			
            // check euc dists btwn FSRs are >= minSelfDistEuc
            SLsquared = Pow(b, 2) + Pow(h, 2); 
			if (SLsquared < Pow(minSL, 2)) { continue; }

            R1 = MakePoint(0.0, 0.0); 
            R2 = MakePoint(l, 0.0); 
            R3 = MakePoint(l - b, h); 
            R4 = MakePoint(-b, h); 

            alpha = atan(h.value / b.value); // bad rigor 

            // choose FSG to lie on boundary of R4 
            for (theta.value = 0; theta < alpha; theta = theta + 2.0 * width) {
                counter++; 
                FSG.x = -b + cos(theta.value); 
                FSG.y = h - sin(theta.value);
                // printf("FSG = (%6.5f, %6.5f). \n", FSG.x.value, FSG.y.value); 

                // check that FSG is not too close to R3 and R1
                if (distSquared(FSG, R3) < 1.0 || distSquared(FSG, R1) < 1.0) { continue; }

                intersectR1 = isIntersectTwoCirc(FSG, R1, distGG.value, distRG.value);
                intersectR2 = isIntersectTwoCirc(FSG, R2, distGG.value, distRG.value);  
                intersectR3 = isIntersectTwoCirc(FSG, R3, distGG.value, distRG.value);  
                // intersectR4 = isIntersectTwoCirc(FSG, R4, distGG.value, distRG.value);  // should always be true for FSG lying on R4

                // compute arcs of Circ(FSG, distGG) which lie outside of disks centered at the Ri with radius distRG, working clockwise from R4
                // for each arc:
                    // determine whether it is covered by the involved Ri, R(i + 1)
                    // determine whether it is covered by translates of Circ(FSG, distGG), working clockwise

                temp1 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value); 

                if (intersectR3) {
                    // ("    checking arc R4 to R3.\n"); 
                    temp2 = intersectTwoCircGen1(FSG, R3, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = acos((temp2.x.value - FSG.x.value)/ distGG.value);

                    // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, do nothing
                    if (alpha1 > alpha2) {
                        // else consider placements of G on arc and verify that 
                        // each placement is covered by disks centered at R1, R2, or translates of FSG
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0 * width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R1, R2, and translates of FSG
                            // print first counterexample found for each l, b, h
                            if (distSquared(G, R1) > distRGsq
                                && distSquared(G, R2) > distRGsq
                                && distSquared(G, FSG) > distGGsq
                                && distSquared(G, FSG - R4) > distGGsq
                                && distSquared(G, FSG - R2) > distGGsq
                                && distSquared(G, FSG + R4) > distGGsq
                                && distSquared(G, FSG + R2) > distGGsq) {
                                
                                printf("1, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            } 

                            else {
                                // printf("no configs so far.\n"); 
                            }
                        }
                        
                    }
                    
                    else {
                        // printf("1 no configs so far.\n"); 
                    }
                }

                // The rest of the code is exceptionally repetitive.

                // if Circ(FSG, distGG) does not intersect R3 but intersects R2, check the arc going from R4 to R2
                else if (intersectR2) {
                    // printf("    checking arc R4 to R2.\n");
                    temp2 = intersectTwoCircGen1(FSG, R2, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = asin((temp2.y.value - FSG.y.value)/ distGG.value);

                    // if arc is covered by R4 and R2, continue
                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0 * width) {
                            G.x = FSG.x + cos(phi.value); 
                            G.y = FSG.y + sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R1, R3, and translates of FSG
                            if (distSquared(G, R1) > distRGsq
                                && distSquared(G, R3) > distRGsq
                                && distSquared(G, FSG) > distGGsq
                                && distSquared(G, FSG - R4) > distGGsq
                                && distSquared(G, FSG - R2) > distGGsq
                                && distSquared(G, FSG + R4) > distGGsq
                                && distSquared(G, FSG + R2) > distGGsq) { 
                                
                                printf("2, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            } 
                        }
                    }
                }

                // if Circ(FSG, distGG) does not intersect R3 or R2, but intersects R1, check the arc going from R4 to R1
                // this case probably never happens in real life
                else if (intersectR1) {
                    // printf("    checking arc R4 to R1.\n");
                    temp2 = intersectTwoCircGen1(FSG, R1, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);

                    // if arc is covered by R4 and R1, continue
                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0 * width) {
                            G.x = FSG.x + cos(phi.value); 
                            G.y = FSG.y + sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R2, R3, and translates of FSG
                            if (distSquared(G, R2) > distRGsq
                                && distSquared(G, R3) > distRGsq
                                && distSquared(G, FSG) > distGGsq
                                && distSquared(G, FSG - R4) > distGGsq
                                && distSquared(G, FSG - R2) > distGGsq
                                && distSquared(G, FSG + R4) > distGGsq
                                && distSquared(G, FSG + R2) > distGGsq) { 
                                
                                printf("3, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            } 
                        }
                    }

                }

                // also should add the case where Circ(FSG, distGG) does not intersect R3, R2, or R1, but this probably never happens
                else {
                    // printf("    checking arc R4 to R4.\n");
                    temp2 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);
                    
                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200 *width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R1, R2, R3, and translates of FSG
                            if (distSquared(G, R1) > distRGsq
                                    && distSquared(G, R2) > distRGsq
                                    && distSquared(G, R3) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 
                                
                                printf("4, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  
                        }
                    }
                }


                // next, if Circ(FSG, distGG) intersects R3, we want to deal with arcs going from R3 to R2, R3 to R1, R3 to R4 clockwise
                if (intersectR3) {
                    temp1 = intersectTwoCircGen2(FSG, R3, distGG.value, distRG.value);

                    if (intersectR2) {
                        // printf("    checking arc R3 to R2.\n");
                        temp2 = intersectTwoCircGen1(FSG, R2, distGG.value, distRG.value); 
                        alpha1 = asin((temp1.y.value - FSG.y.value)/ distGG.value);
                        alpha2 = asin((temp2.y.value - FSG.y.value)/ distGG.value);

                        // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, continue
                        if (alpha1 > alpha2) { 
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0 * width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R1, R4, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                    && distSquared(G, R4) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 
                                    
                                    printf("5, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                } 
                            }
                        }
                    }

                    else if (intersectR1) {
                        // printf("    checking arc R3 to R1.\n");
                        temp2 = intersectTwoCircGen1(FSG, R1, distGG.value, distRG.value); 
                        alpha1 = -acos((temp1.x.value - FSG.x.value)/ distGG.value);
                        alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);

                        // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, continue
                        if (alpha1 > alpha2) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R4, and translates of FSG
                                if (distSquared(G, R2) > distRGsq
                                    && distSquared(G, R4) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                    printf("6, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                } 
                            }
                        }
                    }
                    
                    // this case probably doesn't happen in real life either
                    // arc going from R3 to R4
                    else {
                        // printf("    checking arc R3 to R4.\n");
                        temp2 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value); 
                        alpha1 = -acos((temp1.x.value - FSG.x.value)/ distGG.value);
                        alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);

                        // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, continue
                        if (alpha1 > alpha2) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0 * width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R1, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                    && distSquared(G, R2) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                    printf("7, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                } 
                            }
                        }

                    }
                }

                // next, if Circ(FSG, distGG) intersects R2, we want to deal with arcs going from R2 to R1, R2 to R4 clockwise
                if (intersectR2) {
                    temp1 = intersectTwoCircGen2(FSG, R2, distGG.value, distRG.value); 
                    alpha1 = -acos((temp1.x.value - FSG.x.value)/ distGG.value);

                    if (intersectR1) {
                        // printf("    checking arc R2 to R1.\n");
                        temp2 = intersectTwoCircGen2(FSG, R1, distGG.value, distRG.value); 
                        alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);

                        // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, continue
                        if (alpha1 > alpha2) { 
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0 * width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R4, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                    && distSquared(G, R2) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 
                                        
                                    printf("8, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                } 
                            }
                        }
                    }

                    // arc from R2 to R4
                    else {
                        // printf("    checking arc R2 to R4.\n");
                        temp2 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value); 
                        alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);

                        // if arc going from intersection w/ R4 to intersection w/ R3 is covered by R3 & R4, continue
                        if (alpha1 > alpha2) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R1, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                    && distSquared(G, R2) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 
                                        
                                    printf("9, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                } 
                            }
                        }
                    }
                }

                // lastly, if Circ(FSG, distGG) intersect R1, we want to deal with the arc going from R1 to R4 clockwise
                if (intersectR1) {
                    // printf("    checking arc R1 to R4.\n");
                    temp1 = intersectTwoCircGen1(FSG, R1, distGG.value, distRG.value);
                    temp2 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value);
                    alpha1 = -acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = -acos((temp2.x.value - FSG.x.value)/ distGG.value);      

                    if (alpha1 > alpha2) {  
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi  = phi - 200.0*width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R3, R4, and translates of FSG
                            if (distSquared(G, R3) > distRGsq
                                    && distSquared(G, R4) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                printf("10, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  
                        }
                    }
                }
            }

            // printf("finished checking all FSGs lying on R4.\n"); 

            for (theta.value = 0; theta < PI - alpha; theta = theta + 2.0 * width ) {
                FSG.x = ZERO + cos(theta.value); 
                FSG.y = ZERO + sin(theta.value); 

                // check that FSG is "close" to border of valid region
                if (distSquared(FSG, R2) < 0.95 || distSquared(FSG, R4) < 0.95) { continue; }

                // intersectR1 = isIntersectTwoCirc(FSG, R1, distGG.value, distRG.value); // should always be true for FSG lying on R1
                intersectR2 = isIntersectTwoCirc(FSG, R2, distGG.value, distRG.value);  
                intersectR3 = isIntersectTwoCirc(FSG, R3, distGG.value, distRG.value);  
                intersectR4 = isIntersectTwoCirc(FSG, R4, distGG.value, distRG.value);  

                // compute arcs of Circ(FSG, distGG) which lie outside of disks centered at the Ri with radius distRG, working clockwise from R1
                // for each arc:
                    // determine whether it is covered by the involved Ri, R(i + 1)
                    // determine whether it is covered by translates of Circ(FSG, distGG), working clockwise
                
                temp1 = intersectTwoCircGen1(FSG, R1, distGG.value, distRG.value);

                // check arc from R1 to R4
                if (intersectR4) {
                    // printf("    checking arc R1 to R4.\n");

                    temp2 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value);
                    alpha1 = PI.value - asin((temp1.y.value - FSG.y.value)/ distGG.value);
                    alpha2 = PI.value - asin((temp2.y.value - FSG.y.value)/ distGG.value); 

                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R2, R3, and translates of FSG
                            if (distSquared(G, R2) > distRGsq
                                    && distSquared(G, R3) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                printf("11, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  

                        }
                    }
                }

                // else check arc from R1 to R3
                else if (intersectR3) {
                    // printf("    checking arc R1 to R3.\n");
                    temp2 = intersectTwoCircGen1(FSG, R3, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = acos((temp2.x.value - FSG.x.value)/ distGG.value); 

                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R2, R4, and translates of FSG
                            if (distSquared(G, R2) > distRGsq
                                    && distSquared(G, R4) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                printf("12, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  
                        }
                    }
                }

                // else check arc from R1 to R2 
                // this case probably doesn't happen in real life
                else if (intersectR2) {
                    // printf("    checking arc R1 to R2.\n");
                    temp2 = intersectTwoCircGen1(FSG, R2, distGG.value, distRG.value); 
                    alpha1 = acos((temp1.x.value - FSG.x.value)/ distGG.value);
                    alpha2 = acos((temp2.x.value - FSG.x.value)/ distGG.value);

                    if (alpha1 > alpha2) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0 * width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R3, R4, and translates of FSG
                            if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                printf("13, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  
                        }
                    }
                }

                // else check arc from R1 to R1
                // this case probably doesn't happen in real life
                else {
                    // printf("    checking arc R1 to R1.\n");
                    temp2 = intersectTwoCircGen1(FSG, R1, distGG.value, distRG.value); 
                    alpha1 = angle(FSG, temp1); 
                    alpha2 = angle(FSG, temp2);

                    bool isCW = (alpha1 > alpha2); 
                    if (alpha1 < alpha2) {
                        G.x = FSG.x + distGG * cos(alpha1 - width); 
                        G.y = FSG.y + distGG * sin(alpha1 - width); 
                        if (distSquared(G, R1) > distRG * distRG) { 
                            isCW = true; 
                            alpha2 = alpha2 - 2*PI.value; 
                            }
                    }

                    if (isCW) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R2, R3, R4, and translates of FSG
                            if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                printf("14, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  
                        }
                    }

                }


                // next, check arcs starting from R4
                if (intersectR4) {
                    temp1 = intersectTwoCircGen2(FSG, R4, distGG.value, distRG.value); 
                    alpha1 = angle(FSG, temp1); 
                    
                    
                    // check arc from R4 to R3
                    if (intersectR3) {
                        // printf("    checking arc R4 to R3.\n");
                        temp2 = intersectTwoCircGen1(FSG, R3, distGG.value, distRG.value); 
                        alpha2 = angle(FSG, temp2);
                        
                        bool isCW = (alpha1 > alpha2); 
                        if (alpha1 < alpha2) {
                            G.x = FSG.x + distGG * cos(alpha1 - width); 
                            G.y = FSG.y + distGG * sin(alpha1 - width); 
                            if (distSquared(G, R3) > distRG * distRG) { 
                                isCW = true; 
                                alpha2 = alpha2 - 2*PI.value; 
                            }
                        }

                        if (isCW) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R1, R2, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                    printf("15, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                }  
                            }
                        }
                    }

                    // check arc from R4 to R2
                    else if (intersectR2) {
                        // printf("    checking arc R4 to R2.\n");
                        temp2 = intersectTwoCircGen1(FSG, R2, distGG.value, distRG.value); 
                        alpha2 = angle(FSG, temp2); 

                        bool isCW = (alpha1 > alpha2); 
                        if (alpha1 < alpha2) {
                            G.x = FSG.x + distGG * cos(alpha1 - width); 
                            G.y = FSG.y + distGG * sin(alpha1 - width); 
                            if (distSquared(G, R2) > distRG * distRG) { 
                                isCW = true; 
                                alpha2 = alpha2 - 2 * PI.value; 
                            }
                        }
                        
                        if (isCW) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R1, R4, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                    printf("16, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                }  
                            }
                        }
                    }
                }
            
                    // check arc from R4 to R1
                    else {
                         //printf("    checking arc R4 to R1.\n");
                        temp2 = intersectTwoCircGen2(FSG, R1, distGG.value, distRG.value); 
                        alpha2 = angle(FSG, temp2); 

                        bool isCW = (alpha1 > alpha2); 
                        if (alpha1 < alpha2) {
                            G.x = FSG.x + distGG * cos(alpha1 - width); 
                            G.y = FSG.y + distGG * sin(alpha1 - width); 
                            if (distSquared(G, R1) > distRG * distRG) { 
                                isCW = true; 
                                alpha2 = alpha2 - 2* PI.value; 
                            }
                        }
                        
                        if (alpha1 > alpha2) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R3, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                    printf("17, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                }  
                            }
                        }
                    }
                

                // next, check arc starting from R3
                if (intersectR3) {
                    temp1 = intersectTwoCircGen2(FSG, R3, distGG.value, distRG.value);
                    alpha1 = angle(FSG, temp1); 

                    // check arc going from R3 to R2
                    if (intersectR2) {
                        // printf("    checking arc R3 to R2.\n");
                        temp2 = intersectTwoCircGen1(FSG, R2, distGG.value, distRG.value); 
                        alpha2 = angle(FSG, temp2); 

                        bool isCW = (alpha1 > alpha2); 
                        if (alpha1 < alpha2) {
                            G.x = FSG.x + distGG * cos(alpha1 - width); 
                            G.y = FSG.y + distGG * sin(alpha1 - width); 
                            if (distSquared(G, R2) > distRG * distRG) { 
                                isCW = true; 
                                alpha2 = (alpha2 - 2 * PI).value; 
                            }
                        }
                        
                        if (isCW) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R1, R4, and translates of FSG
                                if (distSquared(G, R1) > distRGsq
                                        && distSquared(G, R2) > distRGsq
                                        && distSquared(G, R3) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 
                                    printf("18, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                }  
                            }
                        }

                    }

                    // else check arc from R3 to R1
                    else {
                        // printf("    checking arc R3 to R1.\n");
                        temp2 = intersectTwoCircGen2(FSG, R1, distGG.value, distRG.value); 
                        alpha2 = angle(FSG, temp2); 

                        bool isCW = (alpha1 > alpha2); 
                        if (alpha1 < alpha2) {
                            G.x = FSG.x + distGG * cos(alpha1 - width); 
                            G.y = FSG.y + distGG * sin(alpha1 - width); 
                            if (distSquared(G, R1) > distRG * distRG) { 
                                isCW = true; 
                                alpha2 = alpha2 - 2 * PI.value; 
                            }
                        }
                        
                        if (isCW) {
                            for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                                G.x = FSG.x + distGG * cos(phi.value); 
                                G.y = FSG.y + distGG * sin(phi.value);

                                // if G does not lie in P, translate so it does
                                if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                                // check distances between G and R2, R4, and translates of FSG
                                if (distSquared(G, R2) > distRGsq
                                        && distSquared(G, R4) > distRGsq
                                        && distSquared(G, FSG) > distGGsq
                                        && distSquared(G, FSG - R4) > distGGsq
                                        && distSquared(G, FSG - R2) > distGGsq
                                        && distSquared(G, FSG + R4) > distGGsq
                                        && distSquared(G, FSG + R2) > distGGsq) { 

                                    printf("19, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                    l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                    existsConfig = true; 
                                    continue;
                                }  

                                else {
                                // printf("no configs so far.\n"); 
                                }
                            }
                        }

                    }
                }

                // lastly, check arc starting R2 to R1
                if (intersectR2) {
                    // printf("    checking arc R2 to R1.\n");
                    temp1 = intersectTwoCircGen2(FSG, R2, distGG.value, distRG.value); 
                    alpha1 = angle(FSG, temp1); 

                    temp2 = intersectTwoCircGen2(FSG, R1, distGG.value, distRG.value); 
                    alpha2 = angle(FSG, temp2); 

                    bool isCW = (alpha1 > alpha2); 
                    if (alpha1 < alpha2) {
                        G.x = FSG.x + distGG * cos(alpha1 - width); 
                        G.y = FSG.y + distGG * sin(alpha1 - width); 
                        if (distSquared(G, R1) > distRG * distRG) { 
                            isCW = true; 
                            alpha2 = alpha2 - 2*PI.value; 
                        }
                    }
                    
                    if (isCW) {
                        for (phi.value = alpha1 - 100 * width; phi > alpha2; phi = phi - 200.0*width) {
                            G.x = FSG.x + distGG * cos(phi.value); 
                            G.y = FSG.y + distGG * sin(phi.value);

                            // if G does not lie in P, translate so it does
                            if (isInP(l, b, h, G) == false) { G = translateInP(l, b, h, G); }

                            // check distances between G and R3, R4, and translates of FSG
                            if (distSquared(G, R1) > distRGsq   
                                    && distSquared(G, R2) > distRGsq
                                    && distSquared(G, R3) > distRGsq
                                    && distSquared(G, R4) > distRGsq
                                    && distSquared(G, FSG) > distGGsq
                                    && distSquared(G, FSG - R4) > distGGsq
                                    && distSquared(G, FSG - R2) > distGGsq
                                    && distSquared(G, FSG + R4) > distGGsq
                                    && distSquared(G, FSG + R2) > distGGsq) { 

                                printf("20, l = %6.5f, b = %6.5f, h = %6.5f, FSG = (%4.2f, %4.2f), G = (%4.2f, %4.2f).\n", 
                                l.value, b.value, h.value, FSG.x.value, FSG.y.value, G.x.value, G.y.value);
                                existsConfig = true; 
                                continue;
                            }  

                            else {
                                // printf("no configs so far.\n"); 
                            }
                        }
                    }

                }

            }

            if (existsConfig == false) {
                // printf("    no configurations.\n"); 
                continue; 
            }

            printf("finished checking all FSGs lying on R1.\n"); 
            
        }
    }

    printf("program terminated. \n");
}
