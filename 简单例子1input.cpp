// liustest.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<iostream>

using namespace std;

#define NUM_INPUT 1
#define NUM_HIDDEN 8
#define NUM_OUTPUT 1
#define STEP 1
#define NUM_NEURAL (NUM_INPUT+NUM_HIDDEN+NUM_OUTPUT)
#define G 0.001

struct chromo
{
	float w[(NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT)*(NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT)];
	float new_neural[NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT];
	float old_neural[NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT];
	float a, b, c, d;
};

float hebb(float x,float y,float a,float b,float c,float d)
{
	return (a*x + b*y + c*x*y + d)*G;
}

void initial_chromo(struct chromo ch)
{
	srand(0);
	ch.a = (float)rand() / (RAND_MAX + 1);
	ch.b = (float)rand() / (RAND_MAX + 1);
	ch.c = (float)rand() / (RAND_MAX + 1);
	ch.d = (float)rand() / (RAND_MAX + 1);
	for (int i = 0; i < NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT; i++)
	{
		ch.new_neural[i] = 0.0;
		ch.old_neural[i] = 0.5;
	}
	for (int i = 0; i < NUM_NEURAL*NUM_NEURAL; i++)
	{
		ch.w[i] = (float)rand() / (RAND_MAX + 1)*6-3;
	}
}

int main()
{	
	time_t ti;
	srand((unsigned int)time(&ti));
	FILE * p;
	if ((p = fopen("d:\\test\\p.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	struct chromo ch;
	ch.a = (float)rand() / (RAND_MAX + 1)*2-1;
	ch.b = (float)rand() / (RAND_MAX + 1)*2-1;
	ch.c = (float)rand() / (RAND_MAX + 1)*2-1;
	ch.d = (float)rand() / (RAND_MAX + 1)*2-1;
	for (int i = 0; i < NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT; i++)
	{
		ch.new_neural[i] = 0.0;
		ch.old_neural[i] = 0.5;
	}
	for (int i = 0; i < NUM_NEURAL*NUM_NEURAL; i++)
	{
		ch.w[i] = (float)rand() / (RAND_MAX + 1) * 6 - 3;
	}
	float temp = 0;
	float inputs[NUM_INPUT];
	float sumwx = 0;
	float po = 0.0;
	for (int path = 0; path < 2000; path++)
	{
		inputs[0] = path / 5000 + temp;
		fprintf(p, "%f\n", po);
		for (int s = 0; s < STEP; s++)
		{
			for (int i = 0; i < (NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT); i++)
			{
				sumwx = 0;
				for (int j = 0;  j < (NUM_HIDDEN + NUM_INPUT + NUM_OUTPUT); j++)
				{
					sumwx = sumwx + ch.old_neural[j] * ch.w[i*(NUM_NEURAL)+j];
				}
				if (i < NUM_INPUT)
				{
					ch.new_neural[i] = inputs[0];
				}
				else
				{
					if (i == NUM_NEURAL - 1)
					{
						ch.new_neural[i] = sin(sumwx);
					}
					else
					{
						ch.new_neural[i] = exp(-sumwx*sumwx);
					}
					//ch.new_neural[i] = 1 / (1 + exp(-1 * sumwx));
				}
			}

			for (int i = 0; i < NUM_NEURAL; i++)
			{
				for (int j = 0; j < NUM_NEURAL; j++)
				{
					ch.w[i*NUM_NEURAL + j] = ch.w[i*NUM_NEURAL + j] + hebb(ch.old_neural[i], ch.old_neural[j], ch.a, ch.b, ch.c, ch.d);
				}
			}
			for (int i = 0; i < NUM_NEURAL; i++)
			{
				ch.old_neural[i] = ch.new_neural[i];
				ch.new_neural[i] = 0.0;
			}
		}
		po = po + 0.01*ch.old_neural[NUM_NEURAL - 1];
		cout << po << endl;
	}
	fclose(p);
    return 0;
}

