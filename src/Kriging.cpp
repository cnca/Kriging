//Standard Output
#include <stdio.h>
//Standard and Math
#include <stdlib.h>
#include <math.h>
//Matrix Handling
#include "Matrix.h"
#include "Vector3.h"
//File Handling
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
//Time
#include <ctime>

using namespace std;

#define SPHERICAL 0
#define EXPONENTIAL 1
#define GAUSSIAN 2
#define WAVE 3
#define RATIONAL_Q 4
#define CIRCULAR 5

struct Surface {
	vector3 * S;
	int K;
	double minX;
	double minY;
	double maxX;
	double maxY;
};

struct VariogramModel {
	double C0; //NUGGET
	double CX; //SILL (C0+C)
	double A;  //RANGE (Max distance to consider v(a)=SILL)
	int VAR;   //Type of variogram to use
};

std::clock_t start;

//Safely remove from memory dynamically created 2D double array
//a Array to delete
//N first dimension of the array
void delete2DArray(double** a, int N)
{
	for(int i= 0; i > N; ++i)
	{
		delete a[i];
	}
	delete a;
}

//Create in dynamic memory 2D double array
//N first dimension of the array
//M second dimension of the array
double** new2DArray(int N,int M)
{
	printf("create new array N:%d M:%d\n", N, M);
	double ** A = new double*[N];
	for(int i = 0; i < N; ++i)
	{
		A[i] = new double[M];
	}
	return A;
}

/** Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
  * Variogram Model: Spherical.
  * H distance between points to covariate
  * C0 nugget:
  * CX sill: C0+C
  * A range:
  */
double sphericalVariogram(double H, double C0, double CX, double A)
{
	if(H > A)
	{
		return CX;
	}
	else
	{
		return C0 + ((CX-C0) * ((3*H)/(2*A))-(H*H*H)/(2*A*A*A));
	}
}

double exponentialVariogram(double H, double C0, double CX, double A)
{
	return C0 + ((CX-C0) * (1-exp((-H/A))));
}

double gaussianVariogram(double H, double C0, double CX, double A)
{
	double e = (H/A);
	e = e*e;
	return C0 + ((CX-C0) * (1-exp(-e)));
}

double waveVariogram(double H, double C0, double CX, double A)
{
	double e = H/A;
	return C0 + ((CX-C0) * (1-(sin(e)/e)));
}

double rationalQuadraticVariogram(double H, double C0, double CX, double A)
{
	double e = (H*H)/(A*A);
	return C0 + ((CX-C0) * (e/(1+e)));
}

double circularVariogram(double H, double C0, double CX, double A)
{
	if(H > A)
	{
		return CX;
	}
	else
	{
		double e = H/A;
		double p = 2/M_PI;
		double r = sqrt(1-e*e);
		return C0 + ((CX-C0) * (1-p*acos(e)+p*e*r));
	}
}

/** Calculate the distance between the members of a block of xyz points (B) to the members of another block of xyz points (A) based on their planar xy coordinates
  * Store their result in D
  * A, B arrays of points
  *	A/B[i] is the form of {x,y,z}
  *	z might be unknown
  * startA/startB is the first row to calculate, it is assumed that start >= 0
  * endA/endB is the last row to calculate, it is assumed that endA/endB < sizeA/sizeB,
  */
matrix calculateVariogram(vector3 A[], vector3 B[], int M, int N, int startA, int endA, int startB, int endB, VariogramModel v)
{
	matrix D = matrix(M,N);
	for(int a = startA; a < endA; ++a)
	{
		for(int b = startB; b < endB; ++b)
		{
			double d = vector3::distance(A[a],B[b]);
			if(v.VAR == SPHERICAL)
			{
				D.set(a,b,sphericalVariogram(d, v.C0, v.CX, v.A));
			}
			else if(v.VAR == EXPONENTIAL)
			{
				D.set(a,b,exponentialVariogram(d, v.C0, v.CX, v.A));
			}
			else if(v.VAR == GAUSSIAN)
			{
				D.set(a,b,gaussianVariogram(d, v.C0, v.CX, v.A));
			}
			else if(v.VAR == WAVE)
			{
				D.set(a,b,waveVariogram(d, v.C0, v.CX, v.A));
			}
			else if(v.VAR == RATIONAL_Q)
			{
				D.set(a,b,rationalQuadraticVariogram(d, v.C0, v.CX, v.A));
			}
			else if(v.VAR == CIRCULAR)
			{
				D.set(a,b,circularVariogram(d, v.C0, v.CX, v.A));
			}
		}
	}
	return D;
}

matrix calculateVariogram(vector3 A, vector3 B[], int M, matrix & D, VariogramModel v)
{
	for(int b = 0; b < M; ++b)
	{
		double d = vector3::distance(A,B[b]);
		if(v.VAR == SPHERICAL)
		{
			D.set(b,0,sphericalVariogram(d, v.C0, v.CX, v.A));
		}
		else if(v.VAR == EXPONENTIAL)
		{
			D.set(b,0,exponentialVariogram(d, v.C0, v.CX, v.A));
		}
		else if(v.VAR == GAUSSIAN)
		{
			D.set(b,0,gaussianVariogram(d, v.C0, v.CX, v.A));
		}
		else if(v.VAR == WAVE)
		{
			D.set(b,0,waveVariogram(d, v.C0, v.CX, v.A));
		}
		else if(v.VAR == RATIONAL_Q)
		{
			D.set(b,0,rationalQuadraticVariogram(d, v.C0, v.CX, v.A));
		}
		else if(v.VAR == CIRCULAR)
		{
			D.set(b,0,circularVariogram(d, v.C0, v.CX, v.A));
		}
	}
	return D;
}

/**
  *
  */
void printTime(string message)
{
	float lap = float(clock () - start) /  CLOCKS_PER_SEC;
	cout << message << ": " << lap << " Seconds" << endl;
}

/** Calculate the kriging interpolation for a given block of known points given a target resolution
  * Z list of known points
  * 	Z[i] is the form of {x,y,z}
  *	z is known for all instances in Z
  *	the algorithm does not assume that Z is ordered or that it's layout represents a pattern
  *	the algorithm does not assume that Z is subset of the resulting surface
  * N size of Z (qty of known points)
  * RX target horizontal resolution of the interpolated surface
  * RY target vertical resolution of the interpolated surface
  */
Surface kriging(vector3 Z[], int N, double RX, double RY, VariogramModel VAR)
{
	start = clock();
	Surface S;

	//Get BOUNDS of data minX, maxX, minY, maxY
	double minX = 1000000000, maxX = -1000000000, minY = 1000000000, maxY = -1000000000;
	for(int i = 0; i < N; ++i)
	{
		if(Z[i].x < minX) minX = Z[i].x;
		if(Z[i].x > maxX) maxX = Z[i].x;
		if(Z[i].y < minY) minY = Z[i].y;
		if(Z[i].y > maxY) maxY = Z[i].y;
	}
	S.minX = minX;
	S.minY = minY;
	S.maxX = maxX;
	S.maxY = maxY;

	//Calculate target size in coordinates
	double width = maxX - minX;
	double height = maxY - minY;
	//Calculate target size in points given the spatial resolutions RX, RY, width and height of the target set
	//The size of the target set must include the boundaries of the set given by the known points, that's why an extra row and column are added
	int sizeX = (int)(floor(width / RX)) + 1;
	int sizeY = (int)(floor(height / RY)) + 1;

	cout << "SURFACE DIMENSIONS IN UNITS: " << sizeX << " X " << sizeY << endl;

	S.K = sizeX * sizeY;

	//Create target array and calculate its coordinates
	S.S = new vector3[S.K];

	for(int x = 0; x < sizeX; ++x)
	{
		double coordX = minX + x*RX;
		for(int y = 0; y < sizeY; ++y)
		{
			double coordY = minY + y*RY;
			S.S[x*sizeY + y].x = coordX;
			S.S[x*sizeY + y].y = coordY;
		}
	}

	printTime("PREPARATION");

	//Calculate covariance from each known point to each other
	matrix V = calculateVariogram(Z, Z, N, N, 0, N, 0, N, VAR);

	printTime("Z VARIOGRAM");

	V.invert();
	printTime("V INVERTED");

	//Calculate S
	int n = S.K;
	int q = n/100;
	int p = q;
	for(int i = 0; i < S.K; ++i)
	{
		matrix D = matrix(N,1);
		if(q != 0 && i > p)
		{
			p += q;
			int prog = p/q;
			std::cout << "\r" << prog - 1 << "% completed: ";
			std::cout.flush();
		}
		//Calculate Variogram Model for first point
		calculateVariogram(S.S[i], Z, N, D, VAR);
		matrix W = V.I->multiply(D);
		S.S[i].z = 0;
		for(int n = 0; n < N; ++n)
		{
			S.S[i].z += W(n,0) * Z[n].z;
			//printf("W: %f Z: %f Zo: %f S: %f\n", W[n*N+k], Z[n][2], W[n*N+k] * Z[n][2], S.S[k][2]);
		}
	}

	printTime("KRIGING");
	return S;
}

/**
  *
  */
vector3 * parseXYZFile(char * filename, int &c, int &e)
{
 	ifstream file;
	file.open (filename);

	//Count lines in file
	c = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');

	if(c == 0)
	{
		e = -2;
		return 0; // ERROR: THE FILE HAS NO LINES
	}

	file.clear(); // clear bad state after eof
	file.seekg( 0 ); // reset iterator to beginning of file

	//Create output structure to store the xyz coordinates
	vector3 * Z = new vector3[c];
	//read line by line
	int p = 0;
	std::string line;
	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		double x, y, z;
		if (!(iss >> x >> y >> z))
		{
			e = -1;
			delete [] Z;
			return 0; // ERROR: THE FILE IS NOT WELL FORMATTED
		}
		Z[p].x = x; Z[p].y = y; Z[p].z = z;
		++p;
	}

	file.close();
	return Z;
}

/**
  *
  */
void writeToFile(vector3 S[], int N, const char* filename)
{
	ofstream o;
	o.open (filename);
	for(int i = 0; i < N; ++i)
	{
		o << fixed << S[i].x << " " << S[i].y << " " << S[i].z << "\n";
	}
	o.close();
}

int main(int argc, char **argv)
{
	VariogramModel v;
	v.C0 = 0.1;
	v.CX = 8.5;
	v.A = 25;
	v.VAR = GAUSSIAN;

	if(argc < 2)
	{
		printf("USAGE: Kriging xyz_file [nugget sill range variogram]\n");
		return 0;
	}
	if(argc > 2)
	{
		v.C0 = atof(argv[2]);
		if(v.C0 == 0)
		{
			printf("\nNUGGET CANNOT BE 0\n\n");
			return 1;
		}
	}
	if(argc > 3)
	{
		v.CX = atof(argv[3]);
	}
	if(argc > 4)
	{
		v.A = atof(argv[4]);
	}
	if(argc > 5)
	{
		v.VAR = atoi(argv[5]);
	}

	printf("FILENAME: %s\n", argv[1]);

	//Load file from disk given the route of the file in ARGV[1]
	int c = 0;
	int e = 0;
	vector3 * Z = parseXYZFile(argv[1], c, e);
	if(e == -1)
	{
		printf("\nINPUT FILE IS INCORRECT\n\n");
		return 1;
	}
	else if(e == -2)
	{
		printf("\nINPUT FILE HAS NO LINES\n\n");
		return 1;
	}

	cout << c << endl;

	Surface S = kriging(Z, c, 1, 1, v);

	stringstream ss;
	ss << "Kriging-" << argv[1] << "-" << v.C0 << "-" << v.CX << "-" << v.A << ".xyz";
	writeToFile(S.S, S.K, ss.str().c_str());
	printTime("FILE: - " + ss.str());

//	printf("\n");
//
//	for(int s = 0; s < S.K; ++s)
//	{
//		printf("X: %f Y: %f Z: %f\n", S.S[s][0], S.S[s][1], S.S[s][2]);
//	}

	//delete2DArray(S.S,S.K);
	//delete2DArray(Z,c);

	return 0;
}
