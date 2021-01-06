// hebb_nn_on_cpu.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>

#define DIMENSION 1
#define NUM_NEURAL_INPUT 1+DIMENSION
#define NUM_NEURAL_HIDDEN 5
#define NUM_NEURAL_OUTPUT DIMENSION
#define EPOCH 1
#define THREAD_EPOCH 200
#define STEPS 10

#define MUTATION_P  0.03
#define CROSS_P 0.7
#define PER_BENCHMARK_TIMES 3 
#define NUM_BENCHMARK 10
#define NUM_NEURAL_NETWORK (DIMENSION*PER_BENCHMARK_TIMES*NUM_BENCHMARK)  //每一个染色体上进行的神经网络
#define PI 3.1415926535897932384626433832795029
#define E 2.7182818284590452353602874713526625
#define SINGLE_WEIGHT (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1)
#define MAX 100
#define MIN -100
#define BLOCK 8
#define THREAD 100
#define NUM_CHROMO  2//(BLOCK*THREAD/(NUM_BENCHMARK*PER_BENCHMARK_TIMES))
#define G 0.002
#define DELTA_X_MAX 0.01
#define DELTA_X_MIN -0.01
#define INI_W_MAX 0.1
#define INI_W_MIN -0.1

using namespace std;

float inputs[NUM_NEURAL_INPUT];
float inputs_w[NUM_NEURAL_INPUT];
float inputs_div[DIMENSION][NUM_NEURAL_INPUT];

struct NN
{
	float w[SINGLE_WEIGHT];
	float input_w[NUM_NEURAL_INPUT];
	float old_neural[NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT];
	float new_neural[NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT];
	int num_input;
	int num_output;
};

struct CHROMO
{
	NN nn_initial;
	NN nn[PER_BENCHMARK_TIMES*NUM_BENCHMARK];
	float a[NUM_BENCHMARK*PER_BENCHMARK_TIMES], b[NUM_BENCHMARK*PER_BENCHMARK_TIMES], c[NUM_BENCHMARK*PER_BENCHMARK_TIMES], d[NUM_BENCHMARK*PER_BENCHMARK_TIMES];
	float A, B, C, D;
	float short_dis[NUM_BENCHMARK*PER_BENCHMARK_TIMES];
	float initial_position[NUM_BENCHMARK*PER_BENCHMARK_TIMES][DIMENSION];
	float position[NUM_BENCHMARK*PER_BENCHMARK_TIMES][DIMENSION];
};

struct IN_FITNESS
{
	float f[THREAD_EPOCH];
};

struct IN_POSITION
{
	float p[THREAD_EPOCH][DIMENSION];
};

float hebb(float a, float b, float c, float d, float x, float y)
{
	return G * (a*x + b * y + c * x*y + d);
}

float sigmoid(float x)
{
	if (x < -80)
		return 0;
	else if (x > 80)
		return 1;
	else
		return 1.0F / (1.0F + exp(-1 * x));
}

void initial_chromo(struct CHROMO chromo[NUM_CHROMO])
{
	time_t ti;
	srand((unsigned int)time(&ti));

	for (int i = 0; i < NUM_CHROMO; i++)
	{
		for (int k = 0; k < NUM_NEURAL_INPUT; k++)
		{
			chromo[i].nn_initial.input_w[k] = (float)rand() / RAND_MAX *(INI_W_MAX - INI_W_MIN) - INI_W_MAX;
		}

		for (int k = 0; k < SINGLE_WEIGHT; k++)
		{
			chromo[i].nn_initial.w[k] = (float)rand() / RAND_MAX * (INI_W_MAX - INI_W_MIN) - INI_W_MAX;
		}

		chromo[i].A = (float)rand() / RAND_MAX * (INI_W_MAX - INI_W_MIN) - INI_W_MAX;
		chromo[i].B = (float)rand() / RAND_MAX * (INI_W_MAX - INI_W_MIN) - INI_W_MAX;
		chromo[i].C = (float)rand() / RAND_MAX * (INI_W_MAX - INI_W_MIN) - INI_W_MAX;
		chromo[i].D = (float)rand() / RAND_MAX * (INI_W_MAX - INI_W_MIN) - INI_W_MAX;

		for (int num_fun = 0; num_fun < NUM_BENCHMARK*PER_BENCHMARK_TIMES; num_fun++)
		{
			for (int t = 0; t < SINGLE_WEIGHT; t++)
			{
				chromo[i].nn[num_fun].w[t] = chromo[i].nn_initial.w[t];
			}
			for (int t = 0; t < NUM_NEURAL_INPUT; t++)
			{
				chromo[i].nn[num_fun].input_w[t] = chromo[i].nn_initial.input_w[t];
			}
			for (int k = 0; k < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; k++)
			{
				chromo[i].nn[num_fun].old_neural[k] = 0.5;
				chromo[i].nn[num_fun].new_neural[k] = 0.0;
			}
			chromo[i].nn[num_fun].num_input = 2;
			chromo[i].nn[num_fun].num_output = 1;
			for (int j = 0; j < DIMENSION; j++)
			{
				chromo[i].initial_position[num_fun][j] = (float)rand() / RAND_MAX * (MAX - MIN) - MAX;
				chromo[i].position[num_fun][j] = chromo[i].initial_position[num_fun][j];
			}

			chromo[i].a[num_fun] = chromo[i].A;
			chromo[i].b[num_fun] = chromo[i].B;
			chromo[i].c[num_fun] = chromo[i].C;
			chromo[i].d[num_fun] = chromo[i].D;
		}
	}
}

float benchmark(int funnum, float x[][DIMENSION], int num)    //some simple functions
{
	int a = funnum;
	float f, sum1, sum2;
	int i;
	f = 0.0;
	i = 0;
	sum1 = 0.0;
	sum2 = 0.0;

	switch (a)
	{
	case 0:                     //sphere function
	{
		f = 0.0;
		for (i = 0; i < DIMENSION; i++)
		{
			f += x[num][i] * x[num][i];
			//printf("x[%d][%d]=%f, \t", num,i,x[num][i]);
		}
		return(f);
	}

	case 1:						//elliptic function
	{	f = 0.0;
	for (i = 0; i < DIMENSION; i++)
	{
		f += pow(pow(10.0, 6.0), i / (DIMENSION - 1))*x[num][i] * x[num][i];
		//printf("x[%d][%d]=%f, \t", num,i,x[num][i]);
	}
	return(f);
	}
	case 2:						//rastrigin's function
	{	f = 0.0;
	for (i = 0; i < DIMENSION; i++)
	{
		f += (x[num][i] * x[num][i] - 10.0*cos(2.0*PI*x[num][i]) + 10.0);
	}
	return(f);
	}
	case 3:						//ackley function
	{
		for (i = 0; i < DIMENSION; i++)
		{
			sum1 += x[num][i] * x[num][i];
			sum2 += cos(2.0*PI*x[num][i]);
		}
		sum1 = -0.2*sqrt(sum1 / DIMENSION);
		sum2 = sum2 / DIMENSION;
		f = E - 20.0*exp(sum1) - exp(sum2) + 20.0;

		return(f);
	}
	case 4:					//rosenbrock function
	{
		for (i = 0; i < DIMENSION - 1; i++)
		{
			sum1 = x[num][i] * x[num][i] - x[num][i + 1];
			sum2 = x[num][i] - 1;
			f += 100.0*sum1*sum1 + sum2 * sum2;
		}

		return(f);
	}
	case 5:					//bent cigar function
	{
		for (i = 1; i < DIMENSION; i++)
		{
			f += pow(10.0, 6.0)*x[num][i] * x[num][i];
		}
		f += x[num][0] * x[num][0];
		return f;
	}
	case 6:					//zakharov function
	{
		for (i = 0; i < DIMENSION; i++)
		{
			float z = x[num][i];
			sum1 += pow(z, 2);
			sum2 += 0.5*z;
		}
		f = sum1 + pow(sum2, 2) + pow(sum2, 4);

		return(f);
	}
	case 7:					//schwefel function 
	{
		for (i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < i; j++)
			{
				sum1 += x[num][i] * x[num][i];
			}
			sum1 = sum1 * sum1;
			f += sum1;
		}
		return(f);
	}
	case 8:					//griewank's function
	{
		for (i = 0; i < DIMENSION; i++)
		{
			sum1 += x[num][i] * x[num][i];
			sum2 *= cos(x[num][i] / sqrt(1.0 + i));
		}
		f = sum1 / 4000.0 - sum2 + 1.0;
		return(f);
	}
	default:					//discus function
		f = pow(10.0, 6.0)*x[num][0] * x[num][0];
		for (i = 1; i < DIMENSION; i++)
		{
			f += x[num][i] * x[num][i];
		}
		return(f);
	}

}

void  nn_forward_and_update_paramenter2(struct CHROMO ch[], int num, float f, int num_fun3, int num_epoch, float delta_p[])  ////
{
	FILE *f_old_neural, *f_new_neural, *f_delta_w, *f_old_w, *f_new_w;
	if ((f_old_neural = fopen("d:\\test\\f_old_neural.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_new_neural = fopen("d:\\test\\f_new_neural.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_delta_w = fopen("d:\\test\\f_delta_w.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_old_w = fopen("d:\\test\\f_old_w.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_new_w = fopen("d:\\test\\f_new_w.txt", "w")) == NULL) {
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	
	int t = 0;
	for (int i = 0; i < DIMENSION; i++)
	{
		printf("D: %d  position:%f\n", i, ch[num].position[num_fun3][i]);
	}

	for (int j = 0; j < SINGLE_WEIGHT; j++)
	{
		fprintf(f_old_w, "%f\n", ch[num].nn[num_fun3].w[j]);
	}

	//输入
	inputs[0] = f / num_epoch;
	inputs[1] = ch[num].position[num_fun3][0] / MAX;
	inputs[2] = ch[num].position[num_fun3][1] / MAX;

	float delta_w[SINGLE_WEIGHT];
	for (int i = 0; i < SINGLE_WEIGHT; i++)
	{
		delta_w[i] = 0;
	}
	for (int j = 0; j < (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT); j++)
	{
		fprintf(f_old_neural, "%f\n", ch[num].nn[num_fun3].old_neural[j]);
	}
	t = 0;
	//前向传递
	for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
	{
		for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
		{
			if (i != j)
			{
				ch[num].nn[num_fun3].new_neural[i] += (ch[num].nn[num_fun3].w[t] * ch[num].nn[num_fun3].old_neural[j]);
				t++;
			}
		}
		//printf("dimension %d,neural %d:%f\n",k, i, ch[num].nn[num_fun3][k].new_neural[i]);
		if (i < NUM_NEURAL_INPUT)
		{
			ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i] + ch[num].nn[num_fun3].input_w[i] * inputs[i]);
		}
		else
		{
			ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i]);
		}
		fprintf(f_new_neural, "%f\n", ch[num].nn[num_fun3].new_neural[i]);
	}
	//printf("\n");		
	//计算delta_w
	t = 0;
	for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
	{
		for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
		{
			if (i != j)
			{
				delta_w[t] = hebb(ch[num].a[num_fun3], ch[num].b[num_fun3], ch[num].c[num_fun3], ch[num].d[num_fun3], ch[num].nn[num_fun3].old_neural[j], ch[num].nn[num_fun3].old_neural[i]);
				fprintf(f_delta_w, "%f\n", delta_w[t]);
				t++;
			}
		}
		fprintf(f_delta_w, "\n\n");
	}
	//更新神经元
	for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
	{
		ch[num].nn[num_fun3].old_neural[i] = ch[num].nn[num_fun3].new_neural[i];
		ch[num].nn[num_fun3].new_neural[i] = 0;
	}
	//位置更新赋值
	for (int k = 0; k < DIMENSION; k++)
	{
		delta_p[k] = ch[num].nn[num_fun3].old_neural[k + NUM_NEURAL_INPUT + NUM_NEURAL_HIDDEN];
	}
	//更新权重
	t = 0;
	for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
	{
		for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
		{
			if (i != j)
			{
				ch[num].nn[num_fun3].w[t] += delta_w[t];
				fprintf(f_new_w, "%f\n", ch[num].nn[num_fun3].w[t]);
				t++;
			}
		}
	}
	fclose(f_old_neural);
	fclose(f_new_neural);
	fclose(f_new_w);
	fclose(f_old_w);
	fclose(f_delta_w);
}

void  nn_forward_and_update_paramenter1(struct CHROMO ch[], int num, float f, int num_fun3, int num_epoch, float delta_p[])  ////
{
	FILE *f_old_neural, *f_new_neural, *f_delta_w, *f_old_w, *f_new_w;
	if ((f_old_neural = fopen("d:\\test\\f_old_neural.txt", "a")) == NULL)
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_new_neural = fopen("d:\\test\\f_new_neural.txt", "w")) == NULL) 
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_delta_w = fopen("d:\\test\\f_delta_w.txt", "w")) == NULL) 
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_old_w = fopen("d:\\test\\f_old_w.txt", "w")) == NULL) 
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	if ((f_new_w = fopen("d:\\test\\f_new_w.txt", "w")) == NULL) 
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	int t = 0;
	for (int i = 0; i < DIMENSION; i++)
	{
		printf("D: %d  position:%f\n", i, ch[num].position[num_fun3][i]);
	}

	for (int j = 0; j < SINGLE_WEIGHT; j++)
	{	
		fprintf(f_old_w, "%f\n", ch[num].nn[num_fun3].w[j]);
	}

	for (int i = 0; i<DIMENSION; i++)
	{
		inputs_div[i][0] = f / num_epoch;
		inputs_div[i][1] = ch[num].position[num_fun3][i] / MAX;
		//inputs[i][2] = ch[num].short_dis[num_fun3] / MAX;
	}

	float w[DIMENSION][SINGLE_WEIGHT];
	float delta_w[SINGLE_WEIGHT];
	for (int i = 0; i < SINGLE_WEIGHT; i++)
	{
		delta_w[i] = 0;
	}
	for (int j = 0; j < (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT); j++)
	{
		fprintf(f_old_neural, "%f\t", ch[num].nn[num_fun3].old_neural[j]);
	}
	fprintf(f_old_neural, "\n");
	for (int k = 0; k < DIMENSION; k++)						//前向传递
	{
		t = 0;
		for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
			{
				if (i != j)
				{
					ch[num].nn[num_fun3].new_neural[i] += (ch[num].nn[num_fun3].w[t] * ch[num].nn[num_fun3].old_neural[j]);
					t++;
				}
			}
			//printf("dimension %d,neural %d:%f\n",k, i, ch[num].nn[num_fun3][k].new_neural[i]);
			if (i<NUM_NEURAL_INPUT)
			{
				ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i] + ch[num].nn[num_fun3].input_w[i] * inputs_div[k][i]);
			}
			else
			{
				ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i]);
			}
			fprintf(f_new_neural, "%f\n", ch[num].nn[num_fun3].new_neural[i]);
		}

		//printf("\n");
		t = 0;
		fprintf(f_delta_w, "DIMENSION:%d\n\n", k);
		for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
			{
			if (i != j)
			{
				w[k][t] = hebb(ch[num].a[num_fun3], ch[num].b[num_fun3], ch[num].c[num_fun3], ch[num].d[num_fun3], ch[num].nn[num_fun3].old_neural[j], ch[num].nn[num_fun3].old_neural[i]);
				fprintf(f_delta_w, "%f\n", w[k][t]);
				t++;
			}
			}
			fprintf(f_delta_w, "\n\n");
		}
		for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			ch[num].nn[num_fun3].old_neural[i] = ch[num].nn[num_fun3].new_neural[i];
			ch[num].nn[num_fun3].new_neural[i] = 0.0;
		}
		delta_p[k] = ch[num].nn[num_fun3].old_neural[k+NUM_NEURAL_HIDDEN+NUM_NEURAL_INPUT];
	}

	for (int i = 0; i < SINGLE_WEIGHT; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			delta_w[i] += w[j][i];
		}
	}
	t = 0;
	for (int i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
	{
		for (int j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
		{
			if (i != j)
			{
				ch[num].nn[num_fun3].w[t] += delta_w[t];
				fprintf(f_new_w, "%f\n", ch[num].nn[num_fun3].w[t]);
				t++;
			}
		}
	}
	fclose(f_old_neural);
	fclose(f_new_neural);
	fclose(f_new_w);
	fclose(f_old_w);
	fclose(f_delta_w);
}

float thread_sort(int sort_num, struct IN_FITNESS in_f[], int num, float p)
{
	float a = 0, count = 0;
	for (int i = 0; i < sort_num - 1; i++)
	{
		for (int j = 0; j < sort_num - i - 1; j++)
		{
			if (in_f[num].f[j] > in_f[num].f[j + 1])
			{
				a = in_f[num].f[j];
				in_f[num].f[j] = in_f[num].f[j + 1];
				in_f[num].f[j + 1] = a;
			}
		}
	}
	for (int i = 0; i < sort_num; i++)
	{
		if (p == in_f[num].f[i])
		{
			return count + 1;
		}
		else
		{
			count++;
		}
	}
}

int main()
{

	CHROMO ch[2];
	initial_chromo(ch);
	float value_position;
	IN_POSITION in_p[1];
	IN_FITNESS in_f[1];
	int num_epoch = 0;
	int num = 0;
	float ff = 0;
	float in_fitness = 0.0;
	float delta_p[DIMENSION] = { 0 };
	FILE *f_position;
	if ((f_position = fopen("d:\\test\\f_position.txt", "w")) == NULL) 
	{
		printf("Cannot open the file,strike any key to exit!\n");
		getchar();
		exit(0);
	}
	for (num_epoch = 0; num_epoch < THREAD_EPOCH; num_epoch++)
	{
		for (int k = 0; k < DIMENSION; k++)
		{
			in_p[0].p[num_epoch][k] = ch[0].position[0][k];
			fprintf(f_position, "%f\n", ch[0].position[0][k]);
		}
		in_f[0].f[num_epoch] = benchmark(num, ch[0].position, 0);
		ff = thread_sort((num_epoch + 1), in_f, 0, in_f[0].f[num_epoch]);
		for (int ttttt = 0; ttttt < STEPS; ttttt++)
		{
			nn_forward_and_update_paramenter1(ch, 0, ff, 0, (num_epoch + 1), delta_p);
		}
		for (int i = 0; i < DIMENSION; i++)
		{
				ch[0].position[0][i] += (delta_p[i] * (DELTA_X_MAX - DELTA_X_MIN) - DELTA_X_MAX)* MAX;	
		}
		in_fitness = benchmark(num, ch[0].position, 0);
		printf("in_fitness[%d]:  %f", num_epoch, in_fitness);
	}
	for (int k = 0; k < DIMENSION; k++)
	{
		fprintf(f_position, "%f\n", ch[0].position[0][k]);
	}

}