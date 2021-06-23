//
//  Created by W. Wang on 2021/6/22.
//  Copyright © 2021 Wei Wang. All rights reserved.
//

#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <math.h>
#include <fstream>
#include <string.h>
#include <random>
#include <ctime>
#include "nrutil.h" 
#include "tools.hpp"

#define N 50
#define fullN_0 160
#define fullN_1 100
#define PI 3.1415926535

double chi(double theta, double  x, double y, float k=0.01);


class Field
{
private:

public:

	static float thickness;
	float r0;
	static float pixel;
	double v0[2];
	double vv;
	double polarity;
	int subdomain[2];
	double center[2];

	std::default_random_engine generator;

	double **phi;
	static float phi0;
	double **fullphi;

    Field();

	Field(double x0, double y0, float r, double *v);


	static void index2position(int i, int j, double *pos);

	void subdomain_shift(double x, double y);

	void shift();

	void fulllattice();

	void fulllattice2subdomain(double ** h, double ** hs);

	static void central_diff(int N1, int N2, double ** phi, int axis, double ** res);

	static void second_diff(int N1, int N2, double ** phi, int axis, double ** res);

	double area();

	void center_of_mass();

	void update_polarity(float dt, bool chemo=false, float k_theta=0.01);

	void update(float dt=0.5, double **h=0, double **h1_laplace=0, double **c=0, 
			bool multi=false, bool adh=false, bool confine=false, bool chemo=false, 
			float k_theta=0.01, float omega=0.05);

	void simulation(float T, float dt=0.5, float ti=0.);

	void save_data(const char * str);
};

class System : public Field
{
private:
	std::default_random_engine generator;

public:
    double **confinement, **conf;
    Field * sys;
    int num_of_field;

    System();

    void save_data(const char *str);

    void update(float dt=0.5, bool therm=false, float k_theta=0.01, float omega=0.05);

    int simulation(const char *str, float T, float dt=0.5, bool therm=false, float k_theta=0.01, float omega=0.05);

};

#endif /* _FIELD_HPP_ */