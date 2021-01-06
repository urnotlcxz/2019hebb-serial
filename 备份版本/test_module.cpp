// test_module.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define NUM_NEURAL_INPUT 2
#define NUM_NEURAL_HIDDEN 2
#define NUM_NEURAL_OUTPUT 1
#define EPOCH 100
#define NUM_CHROMO 100
#define STEPS 10
#define DIMENSION 2
#define MUTATION_P  0.03
#define CROSS_P 0.7
#define PER_BENCHMARK_TIMES 3    
#define NUM_BENCHMARK 10
#define NUM_NEURAL_NETWORK (DIMENSION*PER_BENCHMARK_TIMES*NUM_BENCHMARK)  //每一个染色体上进行的神经网络
#define MY_THREAD 1024
#define PI 3.1415926535897932384626433832795029
#define E 2.7182818284590452353602874713526625
#define SINGLE_WEIGHT (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1)
#define MAX 100
#define MIN -100
struct NN
{
	float w[(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1)];
	float b[(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1)];
	float old_neural[NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT];
	float new_neural[NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT];
	int num_input;
	int num_output;
};

struct CHROMO
{
	NN nn_initial[DIMENSION] = { 0 };
	NN nn[DIMENSION] = {0};
	float A, B, C, D;
	float initial_position[NUM_BENCHMARK*PER_BENCHMARK_TIMES][DIMENSION];
	float position[NUM_BENCHMARK*PER_BENCHMARK_TIMES][DIMENSION];
};

void Sort(float *arr, int length, int *index) {
	int i, j, tin;
	tin = 0;
	float temp = 0.0;
	for (i = 0; i < length - 1; i++)
	{
		for (j = 0; j < length - i - 1; j++)
		{
			if (arr[j]>arr[j + 1])
			{
				temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
				tin = index[j];
				index[j] = index[j + 1];
				index[j + 1] = tin;
			}
		}
	}
}
//锦标赛选择
void selection_tournment(struct CHROMO ch[NUM_CHROMO], float GA_Fitness[NUM_CHROMO], int g)
{
	int m;
	m = NUM_CHROMO;
	int sel_num, sel_init, sel_final;
	sel_init = floor(NUM_CHROMO*0.15);
	sel_final = floor(NUM_CHROMO*0.666);
	sel_num = (int)(floor)(sel_final - sel_init) / (float)EPOCH*g + sel_init;

	int *select;
	select = (int(*))malloc(sel_num * sizeof(int));
	int mark2[NUM_CHROMO] = { 0 };//标记个体有没有被选中
	int index[NUM_CHROMO] = { 0 };//
	float min;
	int maxindex;
	int i, j, d,l, n, k;

	for (i = 0; i < m; i++)			// m = CHROMOSOME_NUM
	{
		for (j = 0; j < sel_num; j++)
		{
			int r2 = rand() % m + 1;   //1-m之间哪个个体（整数）
			while (mark2[r2 - 1] == 1)
			{
				r2 = rand() % m + 1;
			}
			mark2[r2 - 1] = 1;
			select[j] = r2 - 1;
		}
		min = GA_Fitness[select[0]];
		maxindex = select[0];
		for (k = 1; k < sel_num; k++)
		{
			if (GA_Fitness[select[k]] < min)
			{
				min = GA_Fitness[select[k]];
				maxindex = select[k];
			}
		}
		index[i] = maxindex;
		for (n = 0; n < NUM_CHROMO; n++)
		{
			mark2[n] = 0;
		}

		for (n = 0; n < sel_num; n++)
		{
			select[n] = 0;
		}
		for ( d = 0; d < DIMENSION; d++)
		{
			for (l = 0; l < SINGLE_WEIGHT; l++)
			{
				ch[i].nn_initial[d].w[l] = ch[index[i]].nn_initial[d].w[l];
				ch[i].nn_initial[d].b[l] = ch[index[i]].nn_initial[d].b[l];
			}
		}
		ch[i].A = ch[index[i]].A;
		ch[i].B = ch[index[i]].B;
		ch[i].C = ch[index[i]].C;
		ch[i].D = ch[index[i]].D;	
	}
	free(select);
} 
//交叉
void crossover(struct CHROMO ch[NUM_CHROMO])
{
	//const double a = 0.0;
	//const double b = 1.0;
	int two;
	int one;
	int first = 0;
	float r, r2;
	int point;
	float t;
	int i,d;
	for (two = 0; two < NUM_CHROMO; two++)
	{
		r = rand() /RAND_MAX;
		if (r < CROSS_P)
		{
			++first;
			if (first % 2 == 0)//交叉
			{
				//point = rand() % TOTAL_WEIGHT + 1;  //随机选择交叉点
				for ( d = 0; d < DIMENSION; d++)
				{
					for (i = 0; i < SINGLE_WEIGHT; i++)
					{
						r2 = rand() / RAND_MAX;
						if (r2 < 0.5)
						{
							t = ch[one].nn_initial[d].w[i];
							ch[one].nn_initial[d].w[i] = ch[two].nn_initial[d].w[i];
						    ch[two].nn_initial[d].w[i] = t;
						}
					}
					for (i = 0; i < SINGLE_WEIGHT; i++)
					{
						r2 = rand() / RAND_MAX;
						if (r2 < 0.5)
						{
							t = ch[one].nn_initial[d].b[i];
							ch[one].nn_initial[d].b[i] = ch[two].nn_initial[d].b[i];
							ch[two].nn_initial[d].b[i] = t;
						}
					}
				}
				r2 = rand() / RAND_MAX;
				if (r2 < 0.5)
				{
					ch[one].A = ch[two].A;
				}
				r2 = rand() / RAND_MAX;
				if (r2 < 0.5)
				{
					ch[one].B = ch[two].B;
				}
				r2 = rand() / RAND_MAX;
				if (r2 < 0.5)
				{
					ch[one].C = ch[two].C;
				}
				r2 = rand() / RAND_MAX;
				if (r2 < 0.5)
				{
					ch[one].D = ch[two].D;
				}			
			}
			else
			{
				one = two;
			}
		}
	}
}
//变异
void mutation(struct CHROMO ch[NUM_CHROMO])
{
	float r;
	int i, j,d;
	for (i = 0; i < NUM_CHROMO; i++)
	{
		for (d = 0; d < DIMENSION; d++)
		{
			for (j = 0; j < SINGLE_WEIGHT; j++)
			{
				r = rand() / (RAND_MAX + 1.0);
				if (r < MUTATION_P)
				{
					ch[i].nn_initial[d].w[j] = (float)rand() / RAND_MAX * 2.0 - 1.0;
					ch[i].nn_initial[d].b[j] = (float)rand() / RAND_MAX * 2.0 - 1.0;
				}
			}
		}
		r = rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].A = (float)rand() / RAND_MAX * 2.0 - 1.0;	
		}
		r = rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].B = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
		r = rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].C = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
		r = rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].D = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}

	}
}

void initial_chromo(struct CHROMO chromo[NUM_CHROMO])
{
	time_t tttt;
	srand((unsigned int)time(&tttt));
	int t = 0, i = 0, j = 0, k = 0;
	for (i = 0; i < NUM_CHROMO; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (t = 0; t < (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1);t++)
			{
				chromo[i].nn[j].w[t] = rand() * 2 - 1;
				chromo[i].nn[j].b[t] = rand() * 2 - 1;
				chromo[i].nn_initial[j].b[t] = chromo[i].nn[j].b[t];
				chromo[i].nn_initial[j].w[t] = chromo[i].nn[j].w[t];
			}
			chromo[i].nn[j].num_input = 2;
			chromo[i].nn[j].num_output = 1;
			for (k = 0; k < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; k++)
			{
				chromo[i].nn[j].old_neural[k] = 0.5;
				chromo[i].nn[j].new_neural[k] = 0;
			}
		}
		chromo[i].A = rand() / RAND_MAX * 2 - 1;
		chromo[i].B = rand() / RAND_MAX * 2 - 1;
		chromo[i].C = rand() / RAND_MAX * 2 - 1;
		chromo[i].D = rand() / RAND_MAX * 2 - 1;
		for (j = 0;j<NUM_BENCHMARK*PER_BENCHMARK_TIMES;j++)
		{
			for (k = 0;k<DIMENSION;k++)
			{
				chromo[i].initial_position[j][k] = rand() * (MAX-MIN) - MAX;
				chromo[i].position[j][k] = chromo[i].initial_position[j][k];
			}
		}
	}
}

void initial_chromo(struct CHROMO chromo)
{
	time_t tttt;
	srand((unsigned int)time(&tttt));
	int t = 0, i = 0, j = 0, k = 0;
	
		for (j = 0; j < DIMENSION; j++)
		{

			for (t = 0; t < (NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT)*(NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1);t++)
			{
				chromo.nn[j].w[t] = rand() * 2 - 1;
				chromo.nn[j].b[t] = rand() * 2 - 1;
				chromo.nn_initial[j].b[t] = chromo.nn[j].b[t];
				chromo.nn_initial[j].w[t] = chromo.nn[j].w[t];
			}
			chromo.nn[j].num_input = 2;
			chromo.nn[j].num_output = 1;
			for (k = 0; k < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; k++)
			{
				chromo.nn[j].old_neural[k] = 0.5;
				chromo.nn[j].new_neural[k] = 0;
			}	
		chromo.A = rand() / RAND_MAX  * 2 - 1;
		chromo.B = rand() / RAND_MAX  * 2 - 1;
		chromo.C = rand() / RAND_MAX  * 2 - 1;
		chromo.D = rand() / RAND_MAX  * 2 - 1;
		for (j = 0;j<NUM_BENCHMARK*PER_BENCHMARK_TIMES;j++)
		{
			for (k = 0;k<DIMENSION;k++)
			{
				chromo.initial_position[j][k] = rand() * (MAX-MIN) - MAX;
				chromo.position[j][k] = chromo.initial_position[j][k];
			}
		}
	}
}

float my_function(struct CHROMO ch[NUM_CHROMO],float fitness[NUM_CHROMO*PER_BENCHMARK_TIMES*NUM_BENCHMARK], struct CHROMO GA_best_ch,int g)
{
	int i, j, k, fun_num;
	int temp_nn_index[NUM_CHROMO*PER_BENCHMARK_TIMES];//下标
	float nn_fitness[NUM_BENCHMARK*PER_BENCHMARK_TIMES][NUM_CHROMO*PER_BENCHMARK_TIMES];  //实际通过cec计算得出的适应值
	float nn_fitness_one_f[NUM_CHROMO];//所有染色体在同一个函数上的排序
	int record_sort_index[NUM_CHROMO];   //1个函数下的NUM_CHROMO*PER_BENCHMARK_TIMES个排名名次
	float record_nn_fitness_sort_index[NUM_BENCHMARK*PER_BENCHMARK_TIMES][NUM_CHROMO];  //10个函数下的CHROMOSOME_NUM个NN 名次
	float fitness_ranking[NUM_BENCHMARK+1][NUM_CHROMO];
	float GA_fitness[NUM_CHROMO];
	float min_fit;
	FILE* f_temp_ch = fopen("temp_bestch.txt", "a");
	FILE* f_temp_result = fopen("temp_result.txt", "w+");
	//排100个，分别排3次的话
	for (fun_num = 0; fun_num < NUM_BENCHMARK*PER_BENCHMARK_TIMES; fun_num++)
	{
		for (i = 0; i < NUM_CHROMO; i++)
		{
			nn_fitness[fun_num][i] = fitness[fun_num + i*NUM_BENCHMARK*PER_BENCHMARK_TIMES];
		}
	}
	//10个函数分别排3次
	for (fun_num = 0; fun_num < NUM_BENCHMARK*PER_BENCHMARK_TIMES; fun_num++)
	{
		for (i = 0; i < NUM_CHROMO; i++)
		{
			nn_fitness_one_f[i] = nn_fitness[fun_num][i]; //NUM_CHROMO个染色体在一个函数上
			temp_nn_index[i] = i;
		}
		Sort(nn_fitness_one_f, NUM_CHROMO, temp_nn_index);
		for (i = 0; i < NUM_CHROMO; i++)
		{
			record_sort_index[temp_nn_index[i]] = i+1;
		}
		for (i = 0; i < NUM_CHROMO; i++)
		{
			record_nn_fitness_sort_index[fun_num][i] = (float)record_sort_index[i]; //NUM_CHROMO * NUM_BENCHMARK*PER_BENCHMARK_TIMES个名次
		}
	}
	for (i = 0; i < NUM_CHROMO; i++)
	{
		for (fun_num = 0; fun_num < NUM_BENCHMARK; fun_num++)
		{
			for ( k = 0; k < PER_BENCHMARK_TIMES; k++)
			{
				fitness_ranking[fun_num][i] += record_nn_fitness_sort_index[fun_num+ NUM_BENCHMARK*k][i];
			}
			fitness_ranking[fun_num][i] /= PER_BENCHMARK_TIMES;
		}
	}
	//平均排名
	for (i = 0; i < NUM_CHROMO; i++)
	{
		float sum = 0;
		for (j = 0; j < NUM_BENCHMARK; j++)
		{
			sum += fitness_ranking[j][i];
		}
		fitness_ranking[NUM_BENCHMARK][i] = sum / NUM_BENCHMARK;
		GA_fitness[i] = sum / NUM_BENCHMARK;
		printf("the %dth chromosome in all function's avg ranking is :%f \n", i + 1, GA_fitness[i]);
	}
	min_fit = GA_fitness[0];     //当代  找最小排名及对应的网络
	for (k = 0; k < DIMENSION; k++)
	{
		for (j = 0; j < SINGLE_WEIGHT; j++)
		{
			GA_best_ch.nn_initial[k].w[j] = ch[0].nn_initial[k].w[j];
			GA_best_ch.nn_initial[k].b[j] = ch[0].nn_initial[k].b[j];
		}
	}
	GA_best_ch.A = ch[i].A;
	GA_best_ch.B = ch[i].B;
	GA_best_ch.C = ch[i].C;
	GA_best_ch.D = ch[i].D;
	for (i = 1; i < NUM_CHROMO; i++)
	{
		if (GA_fitness[i] < min_fit)
		{
			min_fit = GA_fitness[i];
			for (k = 0; k < DIMENSION; k++)
			{
				for (j = 0; j < SINGLE_WEIGHT; j++)
				{
					GA_best_ch.nn_initial[k].w[j] = ch[i].nn_initial[k].w[j];
					GA_best_ch.nn_initial[k].b[j] = ch[i].nn_initial[k].b[j];
				}
			}
			GA_best_ch.A = ch[i].A;
			GA_best_ch.B = ch[i].B;
			GA_best_ch.C = ch[i].C;
			GA_best_ch.D = ch[i].D;	
		}
	}
	fprintf(f_temp_ch, "g=%d,the best chromosome:\n", g + 1);
	for (k = 0; k < DIMENSION; k++)
	{
		fprintf(f_temp_ch, "The %d dimension \n", k + 1);
		for (j = 0; j < SINGLE_WEIGHT; j++)
		{
			fprintf(f_temp_ch, "%lf\n", GA_best_ch.nn_initial[k].w[j]);
			fprintf(f_temp_ch, "%lf\n", GA_best_ch.nn_initial[k].b[j]);
		}
	}
	fprintf(f_temp_ch, "A:%lf\n B:%lf\n C:%lf\n D:%lf\n", GA_best_ch.A, GA_best_ch.B, GA_best_ch.C, GA_best_ch.D);//当代最佳chromosome
	printf("the generation %d min ranking is:%lf \n", g + 1, min_fit);
	fprintf(f_temp_result, "g=%d,the best:%lf\n", g + 1, min_fit);//当代最佳适应值

	selection_tournment(ch, GA_fitness, g);
	crossover(ch);
	mutation(ch);

	fclose(f_temp_ch);
	fclose(f_temp_result);
	return min_fit;
}

__device__ float hebb(float a, float b, float c, float d, float x, float y)
{
	return (a*x + b * y + c * x*y + d);
}

__device__ float sigmoid(float x)
{
	if (x < -80)
		return 0;
	else if (x > 80)
		return 1;
	else
		return 1.0F / (1.0F + exp(-1 * x));
}

__device__ void  nn_forward_and_update_paramenter(struct CHROMO ch[], int num, float f, int num_fun3)					////
{
	int i, j, k;
	int t = 0;
	for (k = 0; k < DIMENSION; k++)
	{
		ch[num].nn[k].old_neural[0] = f;
		ch[num].nn[k].old_neural[1] = ch[num].position[num_fun3][k];
	}
	for (k = 0; k < DIMENSION; k++)						//前向传递
	{
		t = 0;
		for (i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			for (j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
			{
				if (i != j)
				{
					ch[num].nn[k].new_neural[i] += ch[num].nn[k].w[t] * ch[num].nn[k].old_neural[j] + ch[num].nn[k].b[t];
					t++;
				}
			}
			ch[num].nn[k].new_neural[i] = sigmoid(ch[num].nn[k].new_neural[i] + ch[num].nn[k].old_neural[i]);
		}
	}
	float delta_w = 0;
	for (k = 0; k < DIMENSION; k++)
	{
		t = 0;
		for (i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			for (j = 0; j < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; j++)
			{
				if (i != j)
				{
					delta_w = hebb(ch[num].A, ch[num].B, ch[num].C, ch[num].D, ch[num].nn[k].old_neural[j], ch[num].nn[k].old_neural[i]);
					ch[num].nn[k].w[t] += delta_w;
					ch[num].nn[k].b[t] += delta_w;
					t++;
				}
			}
		}
	}
	for (k = 0; k < DIMENSION; k++)
	{
		for (i = 0; i < NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT; i++)
		{
			ch[num].nn[k].old_neural[i] = ch[num].nn[k].new_neural[i];
			ch[num].nn[k].new_neural[i] = 0;
		}
	}
}

__device__ float benchmark(int funnum, float **x, int num)    //some simple functions
{
	float f, sum1, sum2;
	int i;
	f = 0.0;
	i = 0;
	sum1 = 0.0;
	sum2 = 0.0;
	switch (funnum)
	{
	case 0:                     //sphere function
		for (i = 0; i < DIMENSION; i++)
		{
			f += x[num][i] * x[num][i];
		}
		return(f);
	case 1:						//elliptic function
		for (i = 0; i < DIMENSION; i++)
		{
			f += pow(10.0, 6.0*i / (DIMENSION - 1))*x[num][i] * x[num][i];
		}
		return(f);
	case 2:						//rastrigin's function
		for (i = 0; i < DIMENSION; i++)
		{
			f += (x[num][i] * x[num][i] - 10.0*cos(2.0*PI*x[num][i]) + 10.0);
		}
		return(f);
	case 3:						//ackley function

		for (i = 0; i < DIMENSION; i++)
		{
			sum1 += x[num][i] * x[num][i];
			sum2 += cos(2.0*PI*x[num][i]);
		}
		sum1 = -0.2*sqrt(sum1 / DIMENSION);
		sum2 = sum2 / DIMENSION;
		f = E - 20.0*exp(sum1) - exp(sum2) + 20.0;
		return(f);
	case 4:					//rosenbrock function
		for (i = 0; i < DIMENSION - 1; i++)
		{
			sum1 = x[num][i] * x[num][i] - x[num][i + 1];
			sum2 = x[num][i] - 1;
			f += 100.0*sum1*sum1 + sum2 * sum2;
		}
		return(f);
	case 5:					//bent cigar function
		for (i = 1; i < DIMENSION; i++)
		{
			f += pow(10.0, 6.0)*x[num][i] * x[num][i];
		}
		f += x[num][0] * x[num][0];
		return(f);
	case 6:					//zakharov function
		for (i = 0; i < DIMENSION; i++)
		{
			float z = x[num][i];
			sum1 += pow(z, 2);
			sum2 += 0.5*z;
		}
		f = sum1 + pow(sum2, 2) + pow(sum2, 4);
		return(f);
	case 7:					//schwefel function 
		for (i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < i; j++)
			{
				sum1 += x[num][i] * x[num][i];
			}
			sum1 = sum1*sum1;
			f += sum1;
		}
		return(f);
	case 8:					//griewank's function
		for (i = 0; i < DIMENSION; i++)
		{
			sum1 += x[num][i] * x[num][i];
			sum2 *= cos(x[num][i] / sqrt(1.0 + i));
		}
		f = sum1 / 4000.0 - sum2 + 1.0;
		return(f);
	default:					//discus function
		f = pow(10.0, 6.0)*x[num][0] * x[num][0];
		for (i = 1; i < DIMENSION; i++)
		{
			f += x[num][i] * x[num][i];
		}
		return(f);
	}
}

__global__ void kernel(struct CHROMO ch[NUM_CHROMO], float fitness[NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO])
{
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	if (idx<NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO)
	{
		int num_epoch = 0;
		float ff = 0;
		int num = (idx % (NUM_BENCHMARK*PER_BENCHMARK_TIMES)) % NUM_BENCHMARK;
		do
		{
			ff = benchmark(num, (float **)(ch[idx / (NUM_BENCHMARK*PER_BENCHMARK_TIMES)].position), idx % (NUM_BENCHMARK*PER_BENCHMARK_TIMES));
			for (int ttt = 0; ttt < STEPS; ttt++)
			{
				nn_forward_and_update_paramenter(ch, idx / (NUM_BENCHMARK*PER_BENCHMARK_TIMES), ff, idx % (NUM_BENCHMARK*PER_BENCHMARK_TIMES));
				for (int i = 0; i < DIMENSION; i++)
				{
					ch[idx / (NUM_BENCHMARK*PER_BENCHMARK_TIMES)].position[idx % (NUM_BENCHMARK*PER_BENCHMARK_TIMES)][i] += ch[idx / (NUM_BENCHMARK*PER_BENCHMARK_TIMES)].nn[i].old_neural[NUM_NEURAL_HIDDEN + NUM_NEURAL_INPUT + NUM_NEURAL_OUTPUT - 1];
				}
				ff = benchmark(num, (float **)(ch[idx / (NUM_BENCHMARK*PER_BENCHMARK_TIMES)].position), idx % (NUM_BENCHMARK*PER_BENCHMARK_TIMES));
			}
			num_epoch++;
		} while (num_epoch<200);
		fitness[idx] = ff;
	}
}

int main()
{
	struct CHROMO *h_ch = (struct CHROMO *)malloc(sizeof(CHROMO)*NUM_CHROMO);			//主机端的染色体数组	
	struct CHROMO *d_ch;																//设备端的
	float *d_fitness;																	//设备端的一维适应值数组
	float *h_fitness = (float *)malloc(sizeof(float)*NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO);//主机端的
	int g,d,j,k;
	struct CHROMO GA_best_ch;															//每代最佳的染色体
	struct CHROMO history_best_ch;														//历史最佳的染色体
	float min_fit;																		//每代最佳的适应值
	float the_min_fit_inallG;															//历史最佳的适应值
	int record_g = 0;
	FILE* f_results = fopen("GA200results.txt", "a");
	FILE* f_history_ch = fopen("history_bestch.txt", "a");
	if (f_results == NULL || f_ch == NULL)
	{
		printf("failed to open f file\n");
		system("pause");
	}
	srand((unsigned)time(NULL));

	for (int i = 0;i<NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO;i++)
	{
		h_fitness[i] = 0;
	}
	initial_chromo(h_ch);
	cudaMalloc((void **)&d_fitness, sizeof(float)*NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO);
	cudaMalloc((struct CHROMO **)&d_ch, sizeof(struct CHROMO) * NUM_CHROMO);
	initial_chromo(GA_best_ch);
	initial_chromo(history_best_ch);

	for ( g = 0; g < EPOCH; g++)
	{
		cudaMemcpy(d_ch, h_ch, sizeof(struct CHROMO) * NUM_CHROMO, cudaMemcpyHostToDevice);
		cudaMemcpy(d_fitness, h_fitness, sizeof(float)*NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyHostToDevice);
		kernel << <1, 2 >> >(d_ch, d_fitness);
		cudaMemcpy(d_fitness, h_fitness, sizeof(float)*NUM_BENCHMARK*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyDeviceToHost);
		cudaMemcpy(d_ch, h_ch, sizeof(struct CHROMO)*NUM_CHROMO, cudaMemcpyDeviceToHost);
		min_fit=my_function(h_ch, h_fitness, GA_best_ch,g);
		if (g == 0)
		{
			the_min_fit_inallG = min_fit;
			for (k = 0; k < DIMENSION; k++)
			{
				fprintf(f_history_ch, "The %d dimension \n", k + 1);
				for (j = 0; j < SINGLE_WEIGHT; j++)
				{
					history_best_ch.nn_initial[k].w[j] = GA_best_ch.nn_initial[k].w[j];
					fprintf(f_history_ch, "%lf\n", history_best_ch.nn_initial[k].w[j]);
				}
				for (j = 0; j < SINGLE_WEIGHT; j++)
				{
					history_best_ch.nn_initial[k].b[j] = GA_best_ch.nn_initial[k].b[j];
					fprintf(f_history_ch, "%lf\n", history_best_ch.nn_initial[k].b[j]);
				}

			}
			history_best_ch.A = GA_best_ch.A;
			history_best_ch.B = GA_best_ch.B;
			history_best_ch.C = GA_best_ch.C;
			history_best_ch.D = GA_best_ch.D;
			fprintf(f_history_ch, "A:%lf\n B:%lf\n C:%lf\n D:%lf\n", history_best_ch.A, history_best_ch.B, history_best_ch.C, history_best_ch.D);
			record_g = g;
		}

		if (min_fit < the_min_fit_inallG)
		{
			the_min_fit_inallG = min_fit;
			for (k = 0; k < DIMENSION; k++)
			{
				fprintf(f_history_ch, "The %d dimension \n", k + 1);
				for (j = 0; j < SINGLE_WEIGHT; j++)
				{
					history_best_ch.nn_initial[k].w[j] = GA_best_ch.nn_initial[k].w[j];
					fprintf(f_history_ch, "%lf\n", history_best_ch.nn_initial[k].w[j]);
				}
				for (j = 0; j < SINGLE_WEIGHT; j++)
				{
					history_best_ch.nn_initial[k].b[j] = GA_best_ch.nn_initial[k].b[j];
					fprintf(f_history_ch, "%lf\n", history_best_ch.nn_initial[k].b[j]);
				}

			}
			history_best_ch.A = GA_best_ch.A;
			history_best_ch.B = GA_best_ch.B;
			history_best_ch.C = GA_best_ch.C;
			history_best_ch.D = GA_best_ch.D;
			fprintf(f_history_ch, "A:%lf\n B:%lf\n C:%lf\n D:%lf\n", history_best_ch.A, history_best_ch.B, history_best_ch.C, history_best_ch.D);
			record_g = g;
		}
		printf("the history min ranking is:%lf \n", the_min_fit_inallG);
		fprintf(f_results, "history best is in g=%d:%lf\n", record_g + 1, the_min_fit_inallG);
	}//GA代
	fprintf(fw, "A:%lf\n B:%lf\n C:%lf\n D:%lf\n", history_best_ch.A, history_best_ch.B, history_best_ch.C, history_best_ch.D);
	fclose(f_results);
	cudaFree(d_fitness);
	cudaFree(d_ch);
	free(h_ch);
	free(h_fitness);
    return 0;
}

