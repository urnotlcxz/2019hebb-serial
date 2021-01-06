// test1.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include <stdio.h>
#include<math.h>
#include<stdlib.h>

#define NUM_INPUT_NEURAL 2
#define NUM_HIDDEN_NEURAL 5
#define NUM_OUTPUT_NEURAL 1
#define SINGLE_WEIGHTS (NUM_INPUT_NEURAL+NUM_HIDDEN_NEURAL+NUM_OUTPUT_NEURAL)*(NUM_INPUT_NEURAL+ NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL - 1)
#define DIMENSION 2
#define G 0.002
#define DELTA_W_MAX 0.003
#define DELTA_W_MIN -0.003
#define DELTA_X_MAX 0.1
#define DELTA_X_MIN -0.1
#define EPOCH 300
#define STEPS 10
#define MAX 100
#define MIN -100

struct NN
{
	float old_neural[NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL];
	float new_neural[NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL];
	float w[SINGLE_WEIGHTS];
	float b[SINGLE_WEIGHTS];
	int input;
	int output;
};

struct CHROMO
{
	NN nn[DIMENSION];
	float position[DIMENSION];
	float A, B, C, D;
};

void initial_chromo(CHROMO &ch)
{
	int i;
	float chr[2 * SINGLE_WEIGHTS + 4];
	FILE *read = fopen("history_bestch.txt","r");
	for ( i = 0; i < 2 * SINGLE_WEIGHTS + 4; i++)
	{
		fscanf(read, "%f", &chr[i]);
	}


}

float hebb(float A, float B, float C, float D, float x, float y)
{
	return G*(A*x + B*y + C*x*y + D);
}

void nn_forward_and_update(CHROMO &ch, float f, int num_epoch)
{
	int t = 0;
	for (int k = 0; k < DIMENSION; k++)
	{
		ch.nn[k].old_neural[0] = f / num_epoch;
		ch.nn[k].old_neural[1] = ch.position[k] / 100.0;
	}
	float delta_w[DIMENSION][SINGLE_WEIGHTS];
	float final_delta_w[SINGLE_WEIGHTS];
	for (int i = 0; i < SINGLE_WEIGHTS; i++)
	{
		final_delta_w[i] = 0.0;
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		t = 0;
		for (int i = 0; i < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL); i++)
		{
			for (int j = 0; j < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL);
				j++)
			{
				if (i != j)
				{
					ch.nn[k].new_neural[i] += (ch.nn[k].old_neural[j] * ch.nn[k].w[t]
						+ ch.nn[k].b[t]);
					t++;
				}
			}
			ch.nn[k].new_neural[i] = (ch.nn[k].new_neural[i] + ch.nn[k].old_neural[i]);
		}
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		t = 0;
		for (int i = 0; i < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL); i++)
		{
			for (int j = 0; j < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL);
				j++)
			{
				if (i != j)
				{
					delta_w[k][t] = hebb(ch.A, ch.B, ch.C, ch.D, ch.nn
						[k].old_neural[j], ch.nn[k].old_neural[i]);
					if (delta_w[k][t] > DELTA_W_MAX)
					{
						delta_w[k][t] = DELTA_W_MAX;
					}
					if (delta_w[k][t] < DELTA_W_MIN)
					{
						delta_w[k][t] = DELTA_W_MIN;
					}
					t++;
				}
			}
		}
	}
	for (int i = 0; i < SINGLE_WEIGHTS; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			final_delta_w[i] += delta_w[j][i];
		}
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		t = 0;
		for (int i = 0; i < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL); i++)
		{
			for (int j = 0; j < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL);
				j++)
			{
				if (i != j)
				{
					ch.nn[k].w[t] += final_delta_w[t];
					ch.nn[k].b[t] += final_delta_w[t];
					t++;
				}
			}
		}
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		for (int i = 0; i < (NUM_INPUT_NEURAL + NUM_HIDDEN_NEURAL + NUM_OUTPUT_NEURAL); i++)
		{
			ch.nn[k].old_neural[i] = ch.nn[k].new_neural[i];
			ch.nn[k].new_neural[i] = 0.0;
		}
	}
}

float test_function(float x[])
{
	float value = 0;
	return value;
}

int sort(int num, float fitnesss[], float f)
{
	float value = 0;
	int count = 0;
	for (int i = 0; i < num - 1; i++)
	{
		for (int j = 0; j < num - i - 1; j++)
		{
			if (fitnesss[j] > fitnesss[j + 1])
			{
				value = fitnesss[j];
				fitnesss[j] = fitnesss[j + 1];
				fitnesss[j + 1] = value;
			}
		}
	}
	for (int i = 0; i < num; i++)
	{
		if (f == fitnesss[i])
		{
			return (count + 1);
		}
		else
		{
			count++;
		}
		return 0;
	}

}

int main()
{
	CHROMO ch;

	initial_chromo(ch);
	float fitness[EPOCH];
	int num_epoch = 0;
	float f = 0.0;
	float delta_position = 0.0;
	for (int i = 0; i < EPOCH; i++)
	{
		fitness[i] = 0.0;
	}
	while (num_epoch < EPOCH)
	{
		fitness[num_epoch] = test_function(ch.position);
		f = sort((num_epoch + 1), fitness, fitness[num_epoch]);
		for (int i = 0; i < STEPS; i++)
		{
			nn_forward_and_update(ch, f, (num_epoch + 1));
		}
		for (int k = 0; k < DIMENSION; k++)
		{
			delta_position = ch.nn[k].old_neural[NUM_HIDDEN_NEURAL + NUM_INPUT_NEURAL+ NUM_OUTPUT_NEURAL - 1] * (DELTA_X_MAX - DELTA_X_MIN) - DELTA_X_MAX;
			ch.position[k] += delta_position;
			if (ch.position[k] > MAX)
			{
				ch.position[k] = MAX;
			}
			if (ch.position[k] < MIN)
			{
				ch.position[k] = MIN;
			}
		}
		f = test_function(ch.position);
		printf("fitnesss[%d]: %f\n", num_epoch, f);
		num_epoch++;
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		printf("final position[%d]:  %f\n", ch.position[k]);
	}

}















