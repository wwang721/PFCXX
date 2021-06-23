#include <iostream>
#include <math.h>
#include <cstring>
#include "tools.hpp"

void procBar(int rate, const char *desc)
{

	int ratep = rate/2;
	char bar[52];		//51 “=” s，and one "\0"
	const char *label = "|/-\\";

	memset(bar, 0, 52*sizeof(char));
	
	int i = 0;
	while(i <= ratep)
	{
		bar[i] = '=';
		i++;
	}
	
	printf("%s:[%-51s][%d%%][%c]\r", desc, bar, rate, label[rate%4]);
	fflush(stdout);	

}

double angle_regulate(double theta)
//regulates angles > -pi and <= pi
{
	int n = floor(theta / (2 * PI));
	theta = theta - n * 2 * PI;
	while(theta <= - PI)
		theta += 2 * PI;
	while(theta > PI)
		theta -= 2 * PI;

	return theta;
}
