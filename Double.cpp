/*
 * Double.cpp
 *
 *  Created on: Mar 11, 2023
 *      Author: Marc Lackenby & Rob Meyerhoff
 */

#include "Double.h"

Double MakeDouble(double x, double y)
	{
		Double temp;

		temp.value = x;
		temp.error = y;

		return(temp);
	}

Double operator +(Double x, Double y)
	{
		Double temp;

		temp.value = x.value + y.value;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + (x.error + y.error));

		return(temp);
	}

Double operator +(double x, Double y)
	{
		Double temp;

		temp.value = x + y.value;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + y.error);

		return(temp);
	}

Double operator +(Double x, double y)
	{
		Double temp;

		temp.value = x.value + y;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + x.error);

		return(temp);
	}

Double operator -(Double x, Double y)
	{
		Double temp;

		temp.value = x.value - y.value;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + (x.error + y.error));

		return(temp);
	}

Double operator -(double x, Double y)
	{
		Double temp;

		temp.value = x - y.value;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + y.error);

		return(temp);
	}

Double operator -(Double x, double y)
	{
		Double temp;

		temp.value = x.value - y;
		temp.error = (1.0 + 4.0 * UU) * (fabs(temp.value) * UU + x.error);

		return(temp);
	}

Double operator - (Double x)
	{
		return(MakeDouble(-x.value,x.error));
	}

Double operator *(Double x, Double y)
	{
		Double temp;

		temp.value = x.value * y.value;

		temp.error = (1.0 + 6.0*UU)*((fabs(temp.value)*UU + (x.error * fabs(y.value))) +
								 ((fabs(x.value) * y.error) + (x.error * y.error)));

		return(temp);
	}

Double operator *(double x, Double y)
	{
		Double temp;

		temp.value = x * y.value;

		temp.error = (1.0 + 4.0*UU)*(fabs(temp.value) * UU + fabs(x) * y.error);

		return(temp);
	}

Double operator *(Double x, double y)
	{
		return(y * x);
	}

Double operator /(Double x, Double y)
	{
		Double temp;

		if ((fabs(y.value)-y.error)*fabs(y.value) <= 0)
		{
			printf("Division by zero.\n");
			printf("y.value=%e  y.error=%e\n",y.value,y.error);
			exit(1);
		}

		temp.value = x.value / y.value;
		temp.error = (1.0+4.0*UU) * ( (1.0+6.0*UU)*((fabs(x.value) * y.error +
									             x.error * fabs(y.value))/
									   ((fabs(y.value)-y.error)*fabs(y.value)))
									 + UU * fabs(temp.value)  );

		return(temp);
	}

Double operator /(Double x, double y)
	{
		Double temp;

		if (y == 0.0) printf("Division by zero.\n");

		temp.value = x.value / y;
		temp.error = (1.0+4.0*UU) * ( (1.0+4.0*UU)*((x.error * fabs(y))/
									   (fabs(y)*fabs(y)))
									 + UU * fabs(temp.value)  );


		return(temp);
	}

Double operator /(double x, Double y)
	{
		Double temp, xx;

		xx = MakeDouble(x, 0.0);
		temp =  xx/y;

		return(temp);
	}

Double abs(Double x)
	{
		return(MakeDouble(fabs(x.value),x.error));
	}

Double sqrt(Double x)
	{
		Double temp;
		double tempp;

		tempp = fabs(x.value) - x.error;

		if ((tempp < 0.0) || (x.value < 0.0))
		{
			printf("Square root of negative number.\n");
			printf("x.value=%e  x.error=%e\n",x.value,x.error);
			exit(1);
		}
 		else
 		{
			temp.value = sqrt(x.value);
			temp.error = (1.0+4.0*UU) * ( (1.0+4.0*UU)*temp.value -
										  (1.0-4.0*UU)*sqrt(tempp));
		}

		return(temp);
	}

Double expo(double x, int j)
	{
		Double mid, temp, xx;
		int i;

		xx = MakeDouble(x, 0.0);
		mid = xx;
		temp = MakeDouble(1.0,0.0) + xx;

		for (i = 2; i < j; i++)
		{
			mid = (mid * xx)/i;
			temp = temp + mid;
		}

		return(temp);
	}

Double Pow(Double x, int j)
	{
		Double mid;
		int i;

		mid = x;

		for (i = 1; i < j; i++)
		{
			mid = mid*x;
		}

		return(mid);
	}


double bigDsize(Double x) {
	if (x.error > 0) {
		return (x.value + x.error);
	}
	else { return x.value - x.error; }
	
}

double smallDsize(Double x) {
	if (x.error > 0) {
		return (x.value - x.error);
	}
	else { return x.value + x.error; }
}

int operator <(Double x, Double y)
	{
		return(bigDsize(x) < smallDsize(y));
	}

int operator <(Double x, int y)
	{
		return(bigDsize(x) < y);
	}

int operator <(Double x, double y)
	{
		return(bigDsize(x) < y);
	}

int operator >(Double x, Double y)
	{
		return(smallDsize(x) > bigDsize(y));
	}

int operator >(Double x, int y)
	{
		return(smallDsize(x) > y);
	}

int operator >(Double x, double y)
	{
		return(smallDsize(x) > y);
	}



void printDouble(Double x)
	{
		printf("(%15.12f,%e)",x.value,x.error);
	}

/ poly + sqrt approx of acos(x) for x \in (-1, 1), range is (0, pi)
// using Puiseux series
Double acosApprox(Double x) {
	// make sure input is valid (i.e. crop interval if necessary)
	bool containsONE = false; // tracks if input interval contains +- 1 (approximation blows up in these cases)

	if (bigDsize(x) > 1) {
		containsONE = true;
		x = MakeDouble(0.5 * (1.0 + smallDsize(x)), 0.5 * (1.0 + smallDsize(x)) - smallDsize(x));
	}
	if (smallDsize(x) < -1) {
		containsONE = true;
		x = MakeDouble(0.5 * (-1.0 + bigDsize(x)), 0.5 * (-1.0 + bigDsize(x)) - bigDsize(x));
	}

	if (x > 0 && !containsONE) {
		Double y = 1 - x; 
		Double ysqrt = sqrt(y);
		Double SQRT2 = sqrt(MakeDouble(2.0, UU));
		Double temp = SQRT2 * ysqrt + y * ysqrt / (6 * SQRT2) + 3* y * y * ysqrt / (80 * SQRT2) + 5 * y * y * y * ysqrt / (448 * SQRT2);
		return temp; 
	}

	else if (x > 0 && containsONE) {
		Double y = MakeDouble(1 - smallDsize(x), UU); 
		Double ysqrt = sqrt(y);
		Double SQRT2 = sqrt(MakeDouble(2.0, UU));
		Double temp = SQRT2 * ysqrt + y * ysqrt / (6 * SQRT2) + 3* y * y * ysqrt / (80 * SQRT2) + 5 * y * y * y * ysqrt / (448 * SQRT2);
		temp = MakeDouble(0.5 * temp.value, 0.25 * temp.value);
		// printf("		acos result = %6.5f +- %6.5f.\n", temp.value, temp.error);
		return temp; 
	}

	else if (x < 0 && !containsONE) {
		Double y = 1 + x; 
		Double ysqrt = sqrt(y);
		Double SQRT2 = sqrt(MakeDouble(2.0, UU)); 
		Double temp = SQRT2 * ysqrt + y * ysqrt / (6 * SQRT2) + 3* y * y * ysqrt / (80 * SQRT2) + 5 * y * y * y * ysqrt / (448 * SQRT2);
		return PI - temp; 
	}
	
	else {
		Double y = MakeDouble(1 + bigDsize(x), UU); 
		// printf("		y = %6.5f +- %6.5f.\n", y.value, y.error);
		Double ysqrt = sqrt(y);
		Double SQRT2 = sqrt(MakeDouble(2.0, UU)); 
		Double temp = SQRT2 * ysqrt + y * ysqrt / (6 * SQRT2) + 3* y * y * ysqrt / (80 * SQRT2) + 5 * y * y * y * ysqrt / (448 * SQRT2);
		temp = MakeDouble(0.5 * temp.value, 0.25 * temp.value);
		// printf("		acos result = %6.5f +- %6.5f.\n", temp.value, temp.error);
		return PI - temp; 
	}
	
	// }
}

// poly approx of asin(x) for x \in (-1, 1) using asin(x) = atan(x/sqrt(1 - x^2)), range is (-pi/2, pi/2)
// WARNING: approximation is worse near +-sqrt(1/2) = 0.707... due to atan approximation, but should still be usable
Double asinApprox(Double x) {
	Double temp = sqrt(1 - x* x); 
	temp = atanApprox(x/temp); 
	return temp; 
/*
	int n = 7; // use first n terms of Taylor expansion
	double fact = (2 * n) * (2 * n - 1)/(4.0 * pow(n, 2) * (2*n + 1)); // compute (2n)!/(2^(2n)(n!)^2) * 1/(2n + 1) for error bound
	Double temp = x; 
	for (int i = n - 1; i > 1; i--) {
		temp = x + x * x * temp * (2.0 * i) * (2.0 * i - 1.0)/(4.0 * pow(i, 2)) * (2.0 * i - 1.0)/(2.0 * i + 1.0); 
		fact = fact * (2.0 * i) * (2.0 * i - 1)/(4.0 * pow(i, 2));
	}
	temp = x + x * x * temp / 6.0; 
	fact = fact / 2.0;
	printf("	errfact = %6.5f.\n", fact);
	temp.error = bigDsize(Pow(x, 2 * n + 1)) * fact; 
	return temp; 
*/
}

// WARNING: bounds are extremely bad near +- 1 (specifically interval (+-0.9, +-1.1)), range is (-pi/2, pi/2)
// also, if error of x is large to begin with, will significantly propagate for |x| > 1
Double atanApprox( Double x ) {
    int n = 9; // use first n terms of power series expansion 
    Double xPow = x; // will store value of x^{}
    Double temp;
    
    if (x < -1) { // when x < -1, use -atan(x) = atan(-x) = pi/2 - acot(-x)
        temp = -PI / 2; 
        for (int i = 0; i < n; i++) {
            temp = temp - pow(-1, i) /( (2*i + 1) * xPow ); 
            xPow = xPow * x * x; 
        } 
        // update error
        temp.error = 1 / ((2*n + 1) * smallDsize(xPow)); 
    }

    else if (x > 1) { // when x > 1, use atan(x) = pi/2 - acot(x)
        temp = PI / 2; 
        for (int i = 0; i < n; i++) {
			// printf("boop %d.\n", i);
            temp = temp - pow(-1, i) /( (2*i + 1) * xPow ); 
            xPow = xPow * x * x; 
			printf("xPow %d = %6.5f\n", i, xPow.value); 
			printf("xPowerr %d = %6.5f\n", i, xPow.error); 
        } 
        // update error
        temp.error = 1 / ((2*n + 1) * smallDsize(xPow)); 
		// printf("smallDsize(xPow) = %6.5f", smallDsize(xPow)); 
    }

    else if (x > -1 && x < 1) { 
        temp = MakeDouble(0.0, UU); 
        for (int i = 0; i < n; i++) {
            temp = temp + pow(-1, i) * xPow / (2*i + 1); 
            xPow = xPow * x * x; 
        }
        // update error
        temp.error = bigDsize(xPow) / (2 * n + 1);
    }
    
    // if interval about x contains -1 compute endpoints of interval 
    else if (smallDsize(x) < -1 && bigDsize(x) > -1) {
        double atanLow = -PI.value/2;
        double atanHigh = 0; 

        for (int i = 0; i < n; i++ ) {
            atanLow = (atanLow - pow(-1, i) /( (2*i + 1) * xPow)).value;
            atanHigh = (atanHigh + pow(-1, i) * xPow / (2*i + 1)).value; 
            xPow = xPow * x * x; 
        }
        temp = MakeDouble(0.5 * (atanLow + atanHigh), 0.5 * (atanLow - atanHigh));
        temp.error = bigDsize(xPow) / (2 * n + 1) + 1 / ((2 * n + 1) * smallDsize(xPow)); 
    }
    

    // if interval about x contains 1 ( (smallDsize(x) < 1 && bigDsize(x) > 1) )
    else {
        double atanLow = 0; 
        double atanHigh = PI.value/2; 
        
        for (int i = 0; i < n; i++ ) {
            atanLow = (atanLow + pow(-1, i) * xPow / (2*i + 1)).value; 
            atanHigh = (atanHigh -pow(-1, i) /( (2*i + 1) * xPow)).value;
            xPow = xPow * x * x; 
        }
        temp = MakeDouble(0.5 * (atanLow + atanHigh), 0.5 * (atanLow - atanHigh));
        temp.error = bigDsize(xPow) / (2 * n + 1) + 1 / ((2 * n + 1) * smallDsize(xPow)); 
    }
    return temp;
}

// shifts x (an input to cosApprox or sinApprox) to lie in interval (-pi, pi) 
Double shift( Double x ) {
	if (x < 0) {
		while (x < -PI) {
			x = x + 2 * PI; 
		}
	}

	else {
		while (x > 2*PI) {
			x = x - 2 * PI; 
		}

	}
	return x; 

}

Double cosApprox( Double x ) {
	// first shift x to lie in interval (-pi, pi)
	x = shift(x); 
	int n = 7; // use first n terms of Taylor expansion about 0.0
    int fact = 2*n * (2*n - 1); // compute (2n)!
    Double temp = MakeDouble(1.0, UU);
    for (int i = n; i > 0; --i) {
        temp = 1 - x * x * temp / (2 * i * (2 * i - 1)); 
        fact = fact * (2 * i) * (2 * i - 1); 
    }
    // update error
    temp.error = bigDsize(Pow(x, 2*n)) / fact;
    return temp; 
}

Double sinApprox( Double x ) {
	x = shift(x); 
    int n = 7; // use first n terms of Taylor expansion about 0.0
    int fact = (2 * n + 1) * (2 * n); // compute (2n + 1)!
    Double temp = x; 
    for (int i = n - 1; i > 0; i--) {
        temp = x - x * x * temp / (2 * i * (2 * i + 1));
        fact = fact * (2 * i) * (2 * i + 1); 
    }
    // update error
    temp.error = bigDsize(Pow(x, 2*n + 1))/ fact; 
    return temp; 
}
