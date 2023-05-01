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


double bigDsize(Double x)
	{
		//The formula is not quite right when x.value + x.error is
		//negative, but for the uses of this paper, it's good enough.

		return((1.0 + 4.0*UU) * (x.value + x.error));
	}

double smallDsize(Double x)
	{
		//The formula is not quite right when x.value - x.error is
		//negative, but for the uses of this paper, it's good enough.

		return((1.0 - 4.0*UU) * (x.value - x.error));
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
