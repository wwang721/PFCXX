#include <iostream>
#include "field.hpp"

double chi(double theta, double  x, double y, float k)
{
	double angle = angle_regulate(theta);

	double preferred_angle;
	if(y > -10)
		preferred_angle = PI/2;
	else
		preferred_angle = atan2(-10-y, -x);

	return k*(angle - preferred_angle);
}

float Field::thickness = 2;
float Field::phi0 = 1.;
float Field::pixel = 1.;

Field::Field(){}

Field::Field(double x0, double y0, float r, double *v)
{
	unsigned seed;
	seed = time(0);
	generator.seed(seed);
	
	r0 = r;
	v0[0] = v[0]; 
	v0[1] = v[1];
	vv = sqrt(v0[0] * v0[0] + v0[1] * v0[1]);
	polarity = atan2(v0[1], v0[0]);

	subdomain[0] = floor(x0/pixel);
	subdomain[1] = floor(y0/pixel);
	center[0] = x0 - subdomain[0] * pixel;
	center[1] = y0 - subdomain[1] * pixel;

	double pos[2], rvec[2];

	phi = dmatrix(0, N, 0, N);
	fullphi = dmatrix(0, fullN_0, 0, fullN_1);
    for(int i=0; i<=fullN_0; i++)
        for(int j=0; j<=fullN_1; j++)
            fullphi[i][j] = 0.;

	for(int i=0; i<=N; i++)
	{
		for(int j=0; j<=N; j++)
		{
			index2position(i, j, pos);
			rvec[0] = pos[0] - center[0];
			rvec[1] = pos[1] - center[1];
			double rr = sqrt(rvec[0] * rvec[0] + rvec[1] * rvec[1]);
			phi[i][j] = (phi0/2) * tanh((r0-rr)/thickness) + (phi0/2);
		}
	}
	fulllattice();
}


void Field::index2position(int i, int j, double *pos)
{
	pos[0] = (j - N/2) * pixel;
	pos[1] = (i - N/2) * pixel;
}

void Field::subdomain_shift(double x, double y)
{
	int ishift, jshift;
	jshift = int(floor(x/pixel));
	ishift = int(floor(y/pixel));
	
	double temp[N+1][N+1];
	for(int i=0; i<=N; i++)
	{
		for(int j=0; j<=N; j++)
		{
			int ii, jj;
			ii = i + ishift;
			jj = j + jshift;

			if(ii > N)
			{
				ii = ii - N - 1;
			}
			else if(ii < 0)
			{
				ii = N + ii + 1;
			}

			if(jj > N)
			{
				jj = jj - N - 1;
			}
			else if(jj < 0)
			{
				jj = N + jj + 1;
			}
			temp[i][j] = phi[ii][jj];
		}
	}
	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
			phi[i][j] = temp[i][j];

	subdomain[0] += jshift;
	subdomain[1] += ishift;

	if(subdomain[0] > fullN_1/2)
		subdomain[0] = subdomain[0] - fullN_1;
	else if(subdomain[0] < -fullN_1/2)
		subdomain[0] = subdomain[0] + fullN_1;

	if(subdomain[1] > fullN_0/2)
		subdomain[1] = subdomain[1] - fullN_0;
	else if(subdomain[1] < -fullN_0/2)
		subdomain[1] = subdomain[1] + fullN_0;
}

void Field::shift()
{	
	center_of_mass();
	double d;
	d = sqrt(center[0] * center[0] + center[1] * center[1]);
	if(d > 2 * pixel)
		subdomain_shift(center[0], center[1]);
}

void Field::fulllattice()
{
	for(int i=0; i<=N; i++)
	{
		for(int j=0; j<=N; j++)
		{
			int ii, jj;
			ii = int(i - N/2 + subdomain[1] + fullN_0/2);
			jj = int(j - N/2 + subdomain[0] + fullN_1/2);

			if(ii > fullN_0)
			{
				ii = ii - fullN_0 - 1;
			}
			else if(ii < 0)
			{
				ii = fullN_0 + ii + 1;
			}

			if(jj > fullN_1)
			{
				jj = jj - fullN_1 - 1;
			}
			else if(jj < 0)
			{
				jj = fullN_1 + jj + 1;
			}

			fullphi[ii][jj] = phi[i][j]; 
		}
	}
}

void Field::fulllattice2subdomain(double ** h, double ** hs)
{
	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
		{
			int ii, jj;
			ii = int(i - N/2 + subdomain[1] + fullN_0/2);
			jj = int(j - N/2 + subdomain[0] + fullN_1/2);

			if(ii > fullN_0)
			{
				ii = ii - fullN_0 - 1;
			}
			else if(ii < 0)
			{
				ii = fullN_0 + ii + 1;
			}

			if(jj > fullN_1)
			{
				jj = jj - fullN_1 - 1;
			}
			else if(jj < 0)
			{
				jj = fullN_1 + jj + 1;
			}

			hs[i][j] = h[ii][jj];
		}
}


void Field::central_diff(int N1, int N2, double ** phi, int axis, double ** res)
	/*periodic boundary condition
    */
    //axis = 0: y-axis
    //axis = 1: x-axis
{
	if(axis == 0)
	{
		for(int i=0; i<=N1; i++)
			for(int j=0; j<=N2; j++)
			{
				if(i+1 >= N1)
					res[i][j] = (phi[0][j] - phi[i-1][j])/(2 * pixel);
				else if(i-1 < 0)
					res[i][j] = (phi[i+1][j] - phi[N1-1][j])/(2 * pixel);
				else
					res[i][j] = (phi[i+1][j] - phi[i-1][j])/(2 * pixel);
			}
	}
	else if(axis == 1)
	{
		for(int i=0; i<=N1; i++)
			for(int j=0; j<=N2; j++)
			{
				if(j+1 >= N2)
					res[i][j] = (phi[i][0] - phi[i][j-1])/(2 * pixel);
				else if(j-1 < 0)
					res[i][j] = (phi[i][j+1] - phi[i][N2-1])/(2 * pixel);
				else
					res[i][j] = (phi[i][j+1] - phi[i][j-1])/(2 * pixel);
			}
	}
}

void Field::second_diff(int N1, int N2, double ** phi, int axis, double ** res)
/*periodic boundary condition
    */
    //axis = 0: y-axis
    //axis = 1: x-axis
{
	if(axis == 0)
	{
		for(int i=0; i<=N1; i++)
			for(int j=0; j<=N2; j++)
			{
				if(i+1 >= N1)
					res[i][j] = (phi[0][j] - 2*phi[i][j] + phi[i-1][j])/(pixel * pixel);
				else if(i-1 < 0)
					res[i][j] = (phi[i+1][j] - 2*phi[i][j] + phi[N1-1][j])/(pixel * pixel);
				else
					res[i][j] = (phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(pixel * pixel);
			}
	}
	else if(axis == 1)
	{
		for(int i=0; i<=N1; i++)
			for(int j=0; j<=N2; j++)
			{
				if(j+1 >= N2)
					res[i][j] = (phi[i][0] - 2*phi[i][j] + phi[i][j-1])/(pixel * pixel);
				else if(j-1 < 0)
					res[i][j] = (phi[i][j+1] - 2*phi[i][j] + phi[i][N2-1])/(pixel * pixel);
				else
					res[i][j] = (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(pixel * pixel);
			}
	}
}

double Field::area()
{
	double A(0.);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			double p;
			p = (phi[i][j] + phi[i+1][j] + phi[i][j+1] + phi[i+1][j+1]) / 4;
			A += (p * p) * (pixel * pixel);
		}
	return A;
}

void Field::center_of_mass()
{
	double x(0.), y(0.), m(0.);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			double pos[2], p;
			index2position(i, j, pos);
			p = (phi[i][j] + phi[i+1][j] + phi[i][j+1] + phi[i+1][j+1]) / 4;
			m += p * pixel * pixel;
			x += p * pos[0] * pixel * pixel;
			y += p * pos[1] * pixel * pixel;
		}
	center[0] = (x/m);
	center[1] = (y/m);
}

void Field::update_polarity(float dt, bool chemo, float k_theta)
{
	float mu = 0;
	float D = 0.001;
	float sigma = sqrt(2 * D * dt);
	
	std::normal_distribution<double> normal(mu, sigma);
	polarity += normal(generator);

	if(chemo == true)
	{
		double x, y;
		x = center[0] + subdomain[0]*pixel;
		y = center[1] + subdomain[1]*pixel;
		polarity -= chi(polarity, x, y, k_theta) * dt;
	}
	v0[0] = vv * cos(polarity);
	v0[1] = vv * sin(polarity);
}

void Field::update(float dt, double **h, double **h1_laplace, double **c, 
		bool multi, bool adh, bool confine, bool chemo, float k_theta, float omega)
{
	double **gx, **gy, **laplace, **temp, **hs, **cs, **h1_laplace_sub;
	gx = dmatrix(0, N, 0, N);
	gy = dmatrix(0, N, 0, N);
	central_diff(N, N, phi, 1, gx);
	central_diff(N, N, phi, 0, gy);

	temp = dmatrix(0, N, 0, N);
	laplace = dmatrix(0, N, 0, N);
	second_diff(N, N, phi, 1, temp);
	second_diff(N, N, phi, 0, laplace);
	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
			laplace[i][j] += temp[i][j];

	double deltaV;
	deltaV = 1. - (area() / (PI * ((r0 * phi0)*(r0 * phi0))));

	if(multi == true)
	{
		hs = dmatrix(0, N, 0, N);
		fulllattice2subdomain(h, hs);
	}

	if(confine == true)
	{
		cs = dmatrix(0, N, 0, N);
		fulllattice2subdomain(c, cs);
	}

	if(adh == true)
	{
		h1_laplace_sub = dmatrix(0, N, 0, N);
		fulllattice2subdomain(h1_laplace, h1_laplace_sub);
	}

	float alpha = 1;
	float K = 2;
	float lam = 600;
	float gamma = 10;
	float epsilon = 0.2;

	float epsilon_c = 1.;

	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
		{
			double dphi1, dphi2, dphi3, dphi4, dphi5(0.), dphi6(0.), dphi_adh(0.);

			dphi1 = -(v0[0] * gx[i][j] + v0[1] * gy[i][j]);
			dphi2 = alpha * (phi[i][j]*phi[i][j]*phi[i][j] - 1.5 * (phi[i][j]*phi[i][j]) * phi0
					+ 0.5 * phi[i][j] * phi0 * phi0);
			dphi3 = - K * laplace[i][j];
			dphi4 = - (4 * lam * phi[i][j]/(PI * (r0*phi0) * (r0*phi0))) * deltaV;

			if(multi == true)
				dphi5 = 2 * epsilon * phi[i][j] * (hs[i][j] - phi[i][j] * phi[i][j]);

			if(confine == true)
				dphi6 = 2 * epsilon_c * phi[i][j] * (cs[i][j] * cs[i][j]);

			if(adh == true)
				dphi_adh = -omega * (h1_laplace_sub[i][j] - laplace[i][j]);

			phi[i][j] += (dphi1 - ((dphi2 + dphi3 + dphi4 + dphi5 + dphi6 + dphi_adh)/gamma)) * dt;

		}
    
    free_dmatrix(gx, 0, N, 0, N);
    free_dmatrix(gy, 0, N, 0, N);
    free_dmatrix(laplace, 0, N, 0, N);
    free_dmatrix(temp, 0, N, 0, N);
    free_dmatrix(hs, 0, N, 0, N);
    free_dmatrix(cs, 0, N, 0, N);
    free_dmatrix(h1_laplace_sub, 0, N, 0, N);
    
	shift();
	fulllattice();
	update_polarity(dt, chemo, k_theta);
}

void Field::simulation(float T, float dt, float ti)
{
	float t;
	t = ti;
	int k=0;

	char fileName[50];
    sprintf(fileName, "%s", "data");
	save_data(fileName);

	while(t <= T+ti)
	{
		k += 1;
		t += dt;
		update(dt);
		
		if(k % 20 == 0)
		{
			save_data(fileName);
			procBar(int(k * dt * 100/T));
		}
	}
}

void Field::save_data(const char * str)
{
	char fileName[50];
	sprintf(fileName, "%s.csv", str);

	std::ifstream fin(fileName);
	if(!(fin))
	{
		std::ofstream fout(fileName);
		for(int i=0; i<=fullN_0; i++)
		{
			for(int j=0; j<=fullN_1; j++)
			{
				fout << fullphi[i][j];
			
				if(j!=fullN_1)
				{
					fout << ",";
				}
				else
				{
					fout << "\n";
				}
			}
		}
		fout.close();		
	}
}

/*
int main(int argc, char ** argv)
{
	
	double v[2] = {0.025, 0.025};
	Field a(0, 0, 12, v);
	a.simulation(2000,0.5);
	
	//a.save_data("data");

	return 0;
}
*/

