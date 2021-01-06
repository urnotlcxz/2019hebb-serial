//并行版本输入多维
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define DIMENSION 2
#define INPUT_NEURAL_NUM DIMENSION+1
#define HIDDEN_NEURAL_NUM 5
#define OUTPUT_NEURAL_NUM DIMENSION
#define NEURAL_NUM HIDDEN_NEURAL_NUM + INPUT_NEURAL_NUM + OUTPUT_NEURAL_NUM
#define EPOCH 200
#define THREAD_EPOCH 200
#define STEPS 10

#define MUTATION_P  0.03
#define CROSS_P 0.7
#define PER_BENCHMARK_TIMES 3 
#define BENCHMARK_NUM 10
#define NEURAL_NETWORK_NUM (DIMENSION*PER_BENCHMARK_TIMES*BENCHMARK_NUM)  //每一个染色体上进行的神经网络
#define PI 3.1415926535897932384626433832795029
#define E 2.7182818284590452353602874713526625
#define SINGLE_WEIGHT (NEURAL_NUM)*(NEURAL_NUM - 1)
#define MAX 100
#define MIN -100
#define DEVICE_NUM 1
#define BLOCK 8*4
#define THREAD 192*4
#define NUM_CHROMO  (BLOCK*THREAD/(BENCHMARK_NUM*PER_BENCHMARK_TIMES))
#define G 0.002
#define W_MAX 0.01
#define W_MIN -0.01
#define DELTA_X_MAX 0.1
#define DELTA_X_MIN -0.1


struct NN
{
	float w[SINGLE_WEIGHT];
	float input_w[INPUT_NEURAL_NUM]; ////
	float old_neural[NEURAL_NUM];
	float new_neural[NEURAL_NUM];
};

struct CHROMO
{
	NN nn_initial;
	NN nn[PER_BENCHMARK_TIMES*BENCHMARK_NUM];
	float a[BENCHMARK_NUM*PER_BENCHMARK_TIMES], b[BENCHMARK_NUM*PER_BENCHMARK_TIMES], c[BENCHMARK_NUM*PER_BENCHMARK_TIMES], d[BENCHMARK_NUM*PER_BENCHMARK_TIMES];
	float A, B, C, D;
	float initial_position[BENCHMARK_NUM*PER_BENCHMARK_TIMES][DIMENSION];
	float position[BENCHMARK_NUM*PER_BENCHMARK_TIMES][DIMENSION];
};

struct IN_FITNESS
{
	float f[THREAD_EPOCH];
};

struct IN_POSITION
{
	float p[THREAD_EPOCH][DIMENSION];
};

void Sort(float *arr, int length, int *index)
{
	int i, j, tin;
	tin = 0;
	float temp = 0.0;
	for (i = 0; i < length - 1; i++)
	{
		for (j = 0; j < length - i - 1; j++)
		{
			if (arr[j] > arr[j + 1])
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
	int m, i, j, l, n, k, maxindex;
	m = NUM_CHROMO;
	int sel_num, sel_init, sel_final;
	sel_init = floor(NUM_CHROMO*0.15);
	sel_final = floor(NUM_CHROMO*0.666);
	sel_num = (int)(floor)(sel_final - sel_init) / (float)EPOCH*(g + 1) + sel_init;

	int *select;
	select = (int(*))malloc(sel_num * sizeof(int));
	int mark2[NUM_CHROMO] = { 0 };//标记个体有没有被选中
	int index[NUM_CHROMO] = { 0 };//
	float min;

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
		for (l = 0; l < SINGLE_WEIGHT; l++)
		{
			ch[i].nn_initial.w[l] = ch[index[i]].nn_initial.w[l];
		}
		for (l = 0; l < INPUT_NEURAL_NUM; l++)
		{
			ch[i].nn_initial.input_w[l] = ch[index[i]].nn_initial.input_w[l];
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
	time_t tttt;
	srand((unsigned int)time(&tttt));
	int two;
	int one;
	int first = 0;
	float r, r2;
	//int point;
	float t;
	int i;
	for (two = 0; two < NUM_CHROMO; two++)
	{
		r = (float)rand() / RAND_MAX;
		if (r < CROSS_P)
		{
			++first;
			if (first % 2 == 0)//交叉
			{
				//point = rand() % TOTAL_WEIGHT + 1;  //随机选择交叉点
				for (i = 0; i < SINGLE_WEIGHT; i++)
				{
					r2 = (float)rand() / (RAND_MAX + 1.0);
					if (r2 < 0.5)
					{
						t = ch[one].nn_initial.w[i];
						ch[one].nn_initial.w[i] = ch[two].nn_initial.w[i];
						ch[two].nn_initial.w[i] = t;
					}
				}
				for (i = 0; i < INPUT_NEURAL_NUM; i++)
				{
					r2 = (float)rand() / (RAND_MAX + 1.0);
					if (r2 < 0.5)
					{
						t = ch[one].nn_initial.input_w[i];
						ch[one].nn_initial.input_w[i] = ch[two].nn_initial.input_w[i];
						ch[two].nn_initial.input_w[i] = t;
					}
				}
				r2 = (float)rand() / (RAND_MAX + 1.0);
				if (r2 < 0.5)
				{
					t = ch[one].A;
					ch[one].A = ch[two].A;
					ch[two].A = t;
				}
				r2 = (float)rand() / (RAND_MAX + 1.0);
				if (r2 < 0.5)
				{
					t = ch[one].B;
					ch[one].B = ch[two].B;
					ch[two].B = t;
				}
				r2 = (float)rand() / (RAND_MAX + 1.0);
				if (r2 < 0.5)
				{
					t = ch[one].C;
					ch[one].C = ch[two].C;
					ch[two].C = t;
				}
				r2 = (float)rand() / (RAND_MAX + 1.0);
				if (r2 < 0.5)
				{
					t = ch[one].D;
					ch[one].D = ch[two].D;
					ch[two].D = t;
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
	time_t tttt;
	srand((unsigned int)time(&tttt));
	int i, j;
	for (i = 0; i < NUM_CHROMO; i++)
	{
		for (j = 0; j < INPUT_NEURAL_NUM; j++)
		{
			r = (float)rand() / (RAND_MAX + 1.0);
			if (r < MUTATION_P)
			{
				ch[i].nn_initial.input_w[j] = (float)rand() / (RAND_MAX + 1.0) * (W_MAX - W_MIN) - W_MAX;
			}
		}
		for (j = 0; j < SINGLE_WEIGHT; j++)
		{
			r = (float)rand() / (RAND_MAX + 1.0);
			if (r < MUTATION_P)
			{
				ch[i].nn_initial.w[j] = (float)rand() / (RAND_MAX + 1.0) * (W_MAX - W_MIN) - W_MAX;
			}
		}
		r = (float)rand() / (RAND_MAX + 1.0);
/*
		if (r < MUTATION_P)
		{
			ch[i].A = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
		r = (float)rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].B = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
		r = (float)rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].C = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
		r = (float)rand() / (RAND_MAX + 1.0);
		if (r < MUTATION_P)
		{
			ch[i].D = (float)rand() / RAND_MAX * 2.0 - 1.0;
		}
*/
		if (r < MUTATION_P)
		{
			ch[i].A = (float)rand() / RAND_MAX *  (W_MAX - W_MIN) - W_MAX;
			ch[i].B = (float)rand() / RAND_MAX *  (W_MAX - W_MIN) - W_MAX;
			ch[i].C = (float)rand() / RAND_MAX *  (W_MAX - W_MIN) - W_MAX;
			ch[i].D = (float)rand() / RAND_MAX *  (W_MAX - W_MIN) - W_MAX;
		}

	}
}

void initial_chromo(struct CHROMO chromo[NUM_CHROMO])
{
	time_t ti;
	srand((unsigned int)time(&ti));

	for (int i = 0; i < NUM_CHROMO; i++)
	{
		for (int j = 0; j < INPUT_NEURAL_NUM; j++)
		{
			chromo[i].nn_initial.input_w[j] = (float)rand() / RAND_MAX *(W_MAX - W_MIN) - W_MAX;
		}
		for (int j = 0; j < SINGLE_WEIGHT; j++)
		{
			chromo[i].nn_initial.w[j] = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
		}
		chromo[i].A = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
		chromo[i].B = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
		chromo[i].C = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
		chromo[i].D = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
		for (int num_fun = 0; num_fun < BENCHMARK_NUM*PER_BENCHMARK_TIMES; num_fun++)
		{
			for (int j = 0; j < SINGLE_WEIGHT; j++)
			{
				chromo[i].nn[num_fun].w[j] = chromo[i].nn_initial.w[j];
			}
			for (int j = 0; j < INPUT_NEURAL_NUM; j++)
			{
				chromo[i].nn[num_fun].input_w[j] = chromo[i].nn_initial.input_w[j];
			}
			for (int j = 0; j < NEURAL_NUM; j++)
			{
				chromo[i].nn[num_fun].old_neural[j] = 0.5;
				chromo[i].nn[num_fun].new_neural[j] = 0.0;
			}
			for (int d = 0; d < DIMENSION; d++)
			{
				chromo[i].initial_position[num_fun][d] = (float)rand() / RAND_MAX * (MAX - MIN) - MAX;
				chromo[i].position[num_fun][d] = chromo[i].initial_position[num_fun][d];
			}
			chromo[i].a[num_fun] = chromo[i].A;
			chromo[i].b[num_fun] = chromo[i].B;
			chromo[i].c[num_fun] = chromo[i].C;
			chromo[i].d[num_fun] = chromo[i].D;
		}
	}
}
void initial_chromo(struct CHROMO &chromo)
{
	time_t tttt;
	srand((unsigned int)time(&tttt));
	int i = 0;
	for (i = 0; i<INPUT_NEURAL_NUM; i++)
	{
		chromo.nn_initial.input_w[i] = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	}
	for (i = 0; i<SINGLE_WEIGHT; i++)
	{
		chromo.nn_initial.w[i] = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	}
	chromo.A = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	chromo.B = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	chromo.C = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	chromo.D = (float)rand() / RAND_MAX * (W_MAX - W_MIN) - W_MAX;
	for (int num_fun = 0; num_fun < BENCHMARK_NUM*PER_BENCHMARK_TIMES; num_fun++)
	{
		for (i = 0; i < SINGLE_WEIGHT; i++)
		{
			chromo.nn[num_fun].w[i] = chromo.nn_initial.w[i];
		}
		for (i = 0; i < INPUT_NEURAL_NUM; i++)
		{
			chromo.nn[num_fun].input_w[i] = chromo.nn_initial.input_w[i];
		}
		for (i = 0; i < NEURAL_NUM; i++)
		{
			chromo.nn[num_fun].old_neural[i] = 0.5;
			chromo.nn[num_fun].new_neural[i] = 0.0;
		}
		for (i = 0; i < DIMENSION; i++)
		{
			chromo.initial_position[num_fun][i] = (float)rand() / RAND_MAX * (MAX - MIN) - MAX;
			chromo.position[num_fun][i] = chromo.initial_position[num_fun][i];
		}
		chromo.a[num_fun] = chromo.A;
		chromo.b[num_fun] = chromo.B;
		chromo.c[num_fun] = chromo.C;
		chromo.d[num_fun] = chromo.D;
	}
}

float my_function(struct CHROMO ch[NUM_CHROMO], float fitness[NUM_CHROMO*PER_BENCHMARK_TIMES*BENCHMARK_NUM], struct CHROMO &GA_best_ch, int g)
{
	int temp_nn_index[NUM_CHROMO*PER_BENCHMARK_TIMES];//下标
	float nn_fitness[BENCHMARK_NUM*PER_BENCHMARK_TIMES][NUM_CHROMO];  //实际通过cec计算得出的适应值
	float nn_fitness_one_f[NUM_CHROMO];//所有染色体在同一个函数上的排序
	int record_sort_index[NUM_CHROMO];   //1个函数下的NUM_CHROMO*PER_BENCHMARK_TIMES个排名名次
	float record_nn_fitness_sort_index[BENCHMARK_NUM*PER_BENCHMARK_TIMES][NUM_CHROMO];  //10个函数下的CHROMOSOME_NUM个NN 名次
	float fitness_ranking[BENCHMARK_NUM + 1][NUM_CHROMO];
	float GA_fitness[NUM_CHROMO];
	float min_fit;
	FILE* f_temp_ch;
	FILE* f_temp_result = fopen("temp_result.txt", "a");
	//FILE* all_f_g = fopen("all_f_g.txt", "a");
	FILE* one_g_chromo[NUM_CHROMO + 1];
	char fileName[256];
	for (int i = 0; i < NUM_CHROMO; i++)
	{
		sprintf(fileName, "NN/g%d/ch_%d.txt", g + 1, i);
		one_g_chromo[i] = fopen(fileName, "w");
		if (one_g_chromo[i] == NULL)
		{
			printf("\n Error: Cannot open ch_%d filoe \n",i);
		}
		for (int j = 0; j < INPUT_NEURAL_NUM; j++)
		{
			fprintf(one_g_chromo[i], "%f\n", ch[i].nn_initial.input_w[j]);
		}
		for (int j = 0; j < SINGLE_WEIGHT; j++)
		{
			fprintf(one_g_chromo[i], "%f\n", ch[i].nn_initial.w[j]);
		}
		fprintf(one_g_chromo[i], "A:%f\n B:%f\n C:%f\n D:%f\n", ch[i].A, ch[i].B, ch[i].C, ch[i].D);
		fclose(one_g_chromo[i]);

	}
	for (int i = 0; i < BENCHMARK_NUM + 1; i++)
	{
		for (int j = 0; j < NUM_CHROMO; j++)
		{
			fitness_ranking[i][j] = 0.0;
		}
	}
	for (int fun_num = 0; fun_num < BENCHMARK_NUM*PER_BENCHMARK_TIMES; fun_num++)
	{
		for (int i = 0; i < NUM_CHROMO; i++)
		{
			nn_fitness[fun_num][i] = 0.0;
		}
	}
	//排100个，分别排3次的话
	for (int fun_num = 0; fun_num < BENCHMARK_NUM*PER_BENCHMARK_TIMES; fun_num++)
	{
		for (int i = 0; i < NUM_CHROMO; i++)
		{
			nn_fitness[fun_num][i] = fitness[fun_num + i *BENCHMARK_NUM*PER_BENCHMARK_TIMES];
			//printf("nn_fitness[%d][%d]=%f   \t  ",fun_num,i,nn_fitness[fun_num][i]);
		}
	}
	//10个函数分别排3次
	for (int fun_num = 0; fun_num < BENCHMARK_NUM*PER_BENCHMARK_TIMES; fun_num++)
	{
		for (int i = 0; i < NUM_CHROMO; i++)
		{
			nn_fitness_one_f[i] = nn_fitness[fun_num][i]; //NUM_CHROMO个染色体在一个函数上
			temp_nn_index[i] = i;
		}
		Sort(nn_fitness_one_f, NUM_CHROMO, temp_nn_index);
		for (int i = 0; i < NUM_CHROMO; i++)
		{
			record_sort_index[temp_nn_index[i]] = i + 1;
		}
		for (int i = 0; i < NUM_CHROMO; i++)
		{
			record_nn_fitness_sort_index[fun_num][i] = (float)record_sort_index[i];			//NUM_CHROMO * BENCHMARK_NUM*PER_BENCHMARK_TIMES个名次																							//printf("sort_index[%d][%d]=%f   \t ",fun_num,i,record_nn_fitness_sort_index[fun_num][i]);
		}
	}
	for (int i = 0; i < NUM_CHROMO; i++)
	{
		for (int fun_num = 0; fun_num < BENCHMARK_NUM; fun_num++)
		{
			for (int k = 0; k < PER_BENCHMARK_TIMES; k++)
			{
				fitness_ranking[fun_num][i] += record_nn_fitness_sort_index[fun_num + BENCHMARK_NUM * k][i];
			}
			fitness_ranking[fun_num][i] /= PER_BENCHMARK_TIMES;
		}
	}
	//平均排名
	for (int i = 0; i < NUM_CHROMO; i++)
	{
		float sum = 0;
		for (int j = 0; j < BENCHMARK_NUM; j++)
		{
			sum += fitness_ranking[j][i];
		}
		fitness_ranking[BENCHMARK_NUM][i] = sum / BENCHMARK_NUM;
		GA_fitness[i] = sum / BENCHMARK_NUM;
		printf("the %dth chromosome in all function's avg ranking is :%f \n", i + 1, GA_fitness[i]);
		//fprintf(all_f_g, "g=%d,f[%d]=%f:\n", g + 1, i + 1, GA_fitness[i]);
	}
	sprintf(fileName, "NN/g%d/GA_fitness.txt", g + 1);
	one_g_chromo[NUM_CHROMO] = fopen(fileName, "a");
	if (one_g_chromo[NUM_CHROMO] == NULL)
	{
		printf("\n Error: Cannot open GA_fitness file  \n");
	}
	for (int i = 0; i < NUM_CHROMO; i++)
	{
		fprintf(one_g_chromo[NUM_CHROMO], "%f\n", GA_fitness[i]);
	}
	fclose(one_g_chromo[NUM_CHROMO]);
	min_fit = GA_fitness[0];     //当代  找最小排名及对应的网络
	for (int j = 0; j < INPUT_NEURAL_NUM; j++)
	{
		GA_best_ch.nn_initial.input_w[j] = ch[0].nn_initial.input_w[j];
	}
	for (int j = 0; j < SINGLE_WEIGHT; j++)
	{
		GA_best_ch.nn_initial.w[j] = ch[0].nn_initial.w[j];
	}
	GA_best_ch.A = ch[0].A;
	GA_best_ch.B = ch[0].B;
	GA_best_ch.C = ch[0].C;
	GA_best_ch.D = ch[0].D;
	for (int i = 1; i < NUM_CHROMO; i++)
	{
		if (GA_fitness[i] < min_fit)
		{
			min_fit = GA_fitness[i];
			for (int j = 0; j < INPUT_NEURAL_NUM; j++)
			{
				GA_best_ch.nn_initial.input_w[j] = ch[i].nn_initial.input_w[j];
			}
			for (int j = 0; j < SINGLE_WEIGHT; j++)
			{
				GA_best_ch.nn_initial.w[j] = ch[i].nn_initial.w[j];
			}
			GA_best_ch.A = ch[i].A;
			GA_best_ch.B = ch[i].B;
			GA_best_ch.C = ch[i].C;
			GA_best_ch.D = ch[i].D;
		}
	}
	sprintf(fileName, "NN/g%d/temp_bestch.txt", g + 1);
	f_temp_ch = fopen(fileName, "w");
	fprintf(f_temp_ch, "g=%d,the best chromosome:\n", g + 1);

	for (int j = 0; j < INPUT_NEURAL_NUM; j++)
	{
		fprintf(f_temp_ch, "%f\n", GA_best_ch.nn_initial.input_w[j]);
	}
	for (int j = 0; j < SINGLE_WEIGHT; j++)
	{
		fprintf(f_temp_ch, "%f\n", GA_best_ch.nn_initial.w[j]);
	}
	fprintf(f_temp_ch, "A:%f\n B:%f\n C:%f\n D:%f\n", GA_best_ch.A, GA_best_ch.B, GA_best_ch.C, GA_best_ch.D);//当代最佳chromosome
	printf("the generation %d min ranking is:%f \n", g + 1, min_fit);
	fprintf(f_temp_result, "g=%d,the best:%f\n", g + 1, min_fit);//当代最佳适应值

	selection_tournment(ch, GA_fitness, g);
	crossover(ch);
	mutation(ch);

	for (int i = 0; i < NUM_CHROMO; i++)
	{
		for (int j = 0; j < (PER_BENCHMARK_TIMES*BENCHMARK_NUM); j++)
		{
			for (int d = 0; d < SINGLE_WEIGHT; d++)
			{
				ch[i].nn[j].w[d] = ch[i].nn_initial.w[d];
			}
			for (int k = 0; k < DIMENSION; k++)
			{
				ch[i].position[j][k] = ch[i].initial_position[j][k];
			}
			for (int tt = 0; tt < INPUT_NEURAL_NUM; tt++)
			{
				ch[i].nn[j].input_w[tt] = ch[i].nn_initial.input_w[tt];
			}
			for (int d = 0; d < NEURAL_NUM; d++)
			{
				ch[i].nn[j].old_neural[d] = 0.5;
			}

			ch[i].a[j] = ch[i].A;
			ch[i].b[j] = ch[i].B;
			ch[i].c[j] = ch[i].C;
			ch[i].d[j] = ch[i].D;

		}
	}
	fclose(f_temp_ch);
	fclose(f_temp_result);
	//fclose(all_f_g);
/*
	for (int i = 0; i < NUM_CHROMO + 1; i++)
	{
		fclose(one_g_chromo[i]);
	}
*/
	return min_fit;
}

__device__ float hebb(float a, float b, float c, float d, float x, float y)
{
	return G*(a*x + b * y + c * x*y + d);
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

__device__ void  nn_forward_and_update_paramenter(struct CHROMO ch[], int num, float f, int num_fun3, int num_epoch)  ////
{
	float inputs[INPUT_NEURAL_NUM];
	float delta_w[SINGLE_WEIGHT];
	int t = 0;
	//输入
	inputs[0] = f / num_epoch;
	inputs[1] = ch[num].position[num_fun3][0] / MAX;
	inputs[2] = ch[num].position[num_fun3][1] / MAX;

	for (int i = 0; i < SINGLE_WEIGHT; i++)
	{
		delta_w[i] = 0;
	}
	//前向传递
	for (int i = 0; i < NEURAL_NUM; i++)
	{
		for (int j = 0; j < NEURAL_NUM; j++)
		{
			if (i != j)
			{
				ch[num].nn[num_fun3].new_neural[i] += (ch[num].nn[num_fun3].w[t] * ch[num].nn[num_fun3].old_neural[j]);
				t++;
			}
		}
		//printf("dimension %d,neural %d:%f\n",k, i, ch[num].nn[num_fun3][k].new_neural[i]);
		if (i < INPUT_NEURAL_NUM)
		{
			ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i] + ch[num].nn[num_fun3].input_w[i] * inputs[i]);
		}
		else
		{
			ch[num].nn[num_fun3].new_neural[i] = sigmoid(ch[num].nn[num_fun3].new_neural[i]);
		}
	}
	//printf("\n");		
	//计算delta_w
	t = 0;
	for (int i = 0; i < NEURAL_NUM; i++)
	{
		for (int j = 0; j < NEURAL_NUM; j++)
		{
			if (i != j)
			{
				delta_w[t] = hebb(ch[num].a[num_fun3], ch[num].b[num_fun3], ch[num].c[num_fun3], ch[num].d[num_fun3], ch[num].nn[num_fun3].old_neural[j], ch[num].nn[num_fun3].old_neural[i]);
				t++;
			}
		}
	}
	//更新神经元
	for (int i = 0; i < NEURAL_NUM; i++)
	{
		ch[num].nn[num_fun3].old_neural[i] = ch[num].nn[num_fun3].new_neural[i];
		ch[num].nn[num_fun3].new_neural[i] = 0;
	}
	//更新权重
	t = 0;
	for (int i = 0; i < NEURAL_NUM; i++)
	{
		for (int j = 0; j < NEURAL_NUM; j++)
		{
			if (i != j)
			{
				ch[num].nn[num_fun3].w[t] = (ch[num].nn[num_fun3].w[t] + delta_w[t])*0.99;
				t++;
			}
		}
	}
}

__device__ float benchmark(int funnum, float x[][DIMENSION], int num)    //some simple functions
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
			}
			return(f);
		}
		case 1:						//elliptic function
		{	f = 0.0;
			for (i = 0; i < DIMENSION; i++)
			{
				f += pow(pow(10.0, 6.0), i / (DIMENSION - 1))*x[num][i] * x[num][i];
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

__device__ float thread_sort(int sort_num, struct IN_FITNESS in_f[], int num, float p)
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
	return 0;
}

__global__ void kernel(struct CHROMO ch[NUM_CHROMO], float fitness[BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO], struct IN_FITNESS in_f[BENCHMARK_NUM*NUM_CHROMO*PER_BENCHMARK_TIMES], struct IN_POSITION in_p[BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO])
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int i,d;
	if (idx <BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO)
	{
		float value_position;
		int num_epoch = 0;
		float ff = 0;
		int num = (idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES)) % BENCHMARK_NUM;
		for (num_epoch = 0; num_epoch<THREAD_EPOCH; num_epoch++)
		{
			for (d = 0; d < DIMENSION; d++)
			{
				in_p[idx].p[num_epoch][d] = ch[idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].position[idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES)][d];
			}
			in_f[idx].f[num_epoch] = benchmark(num, ch[idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].position, idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES));
			//printf("ff:%f   \n", in_f[idx].f[num_epoch]);
			ff = thread_sort((num_epoch + 1), in_f, idx, in_f[idx].f[num_epoch]);
			//printf("epoch %d ff: %f\n ", num_epoch, ff);
			for (i = 0; i < STEPS; i++)
			{
				nn_forward_and_update_paramenter(ch, idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES), ff, idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES), (num_epoch + 1));
			}
			for (d = 0; d < DIMENSION; d++)
			{
				ch[idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].position[idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES)][d] += (ch[idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].nn[idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].old_neural[NEURAL_NUM - (d + 1)] * (DELTA_X_MAX - DELTA_X_MIN) - DELTA_X_MAX)* MAX;
			}
		}
		fitness[idx] = benchmark(num, ch[idx / (BENCHMARK_NUM*PER_BENCHMARK_TIMES)].position, idx % (BENCHMARK_NUM*PER_BENCHMARK_TIMES));
		printf("in fitness[%d]: %f\n", idx, fitness[idx]);
	}
}

int main()
{
	struct CHROMO *h_ch = (struct CHROMO *)malloc(sizeof(CHROMO)*NUM_CHROMO);		//主机端的染色体数组
	struct CHROMO *d_ch;					//设备端的
	struct IN_FITNESS *h_in_fitness = (struct IN_FITNESS *)malloc(sizeof(IN_FITNESS)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);
	struct IN_FITNESS *d_in_fitness;
	struct IN_POSITION *h_in_position = (struct IN_POSITION *)malloc(sizeof(IN_POSITION)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);
	struct IN_POSITION *d_in_position = (struct IN_POSITION *)malloc(sizeof(IN_POSITION)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);
	float *d_fitness;					//设备端的一维适应值数组
	float *h_fitness = (float *)malloc(sizeof(float)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);//主机端的
	int i, g, j, d;
	struct CHROMO GA_best_ch;			//每代最佳的染色体
	struct CHROMO history_best_ch;			//历史最佳的染色体
	float min_fit;					//每代最佳的适应值
	float the_min_fit_inallG;			//历史最佳的适应值
	int record_g = 0;
	int device_id = DEVICE_NUM;
	FILE* f_results = fopen("GA200results.txt", "a");
	FILE* f_history_ch = fopen("history_bestch.txt", "a");
	//FILE *f_position = fopen("f_position.txt", "a");
	if (f_results == NULL || f_history_ch == NULL)
	{
		printf("failed to open f file\n");
		system("pause");
	}
	srand((unsigned)time(NULL));

	for (i = 0; i < BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO; i++)
	{
		h_fitness[i] = 0;
	}
	for (i = 0; i < BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO; i++)
	{
		for (j = 0; j < THREAD_EPOCH; j++)
		{
			h_in_fitness[i].f[j] = 0.0;
		}
	}
	for (i = 0; i<BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO; i++)
	{
		for (j = 0; j<THREAD_EPOCH; j++)
		{
			for (d = 0; d<DIMENSION; d++)
			{
				h_in_position[i].p[j][d] = 0.0;
			}
		}
	}
	initial_chromo(h_ch);
	initial_chromo(GA_best_ch);
	initial_chromo(history_best_ch);
	//dim3  grid(BLOCK);
	//dim3  threads(THREAD);
	cudaSetDevice(device_id);
	for (g = 0; g < EPOCH; g++)
	{
		cudaMalloc((void **)&d_fitness, sizeof(float)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);
		cudaMalloc((struct CHROMO **)&d_ch, sizeof(struct CHROMO) * NUM_CHROMO);
		cudaMalloc((struct IN_FITNESS **)&d_in_fitness, sizeof(struct IN_FITNESS)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);
		cudaMalloc((struct IN_POSITION **)&d_in_position, sizeof(struct IN_POSITION)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO);

		cudaMemcpy(d_ch, h_ch, sizeof(struct CHROMO) * NUM_CHROMO, cudaMemcpyHostToDevice);
		cudaMemcpy(d_fitness, h_fitness, sizeof(float)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyHostToDevice);
		cudaMemcpy(d_in_fitness, h_in_fitness, sizeof(struct IN_FITNESS)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyHostToDevice);
		cudaMemcpy(d_in_position, h_in_position, sizeof(struct IN_POSITION)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyHostToDevice);

		kernel <<<BLOCK, THREAD >>> (d_ch, d_fitness, d_in_fitness, d_in_position);
		cudaDeviceSynchronize();

		cudaMemcpy(h_fitness, d_fitness, sizeof(float)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyDeviceToHost);
		cudaMemcpy(h_in_position, d_in_position, sizeof(struct IN_POSITION)*BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO, cudaMemcpyDeviceToHost);
		//cudaMemcpy(h_ch, d_ch, sizeof(struct CHROMO)*NUM_CHROMO, cudaMemcpyDeviceToHost);
		/*
		fprintf(f_position, "第%d代:\n", g + 1);
		for (i = 0; i < BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO; i++)
		{
			fprintf(f_position, "第%d个染色体\n", i + 1);
			for (j = 0; j < THREAD_EPOCH; j++)
			{
				fprintf(f_position, "第%d次移动\t\t\t", j + 1);
				for (d = 0; d < DIMENSION; d++)
				{
					fprintf(f_position, "%f\t\t", h_in_position[i].p[j][d]);
				}
			}
		}
		*/
		/*
		for (i = 0; i < BENCHMARK_NUM*PER_BENCHMARK_TIMES*NUM_CHROMO; i++)
		{
			printf("out fitness[%d] = %f\n", i, h_fitness[i]);
		}
		*/
		min_fit = my_function(h_ch, h_fitness, GA_best_ch, g);
		if (g == 0)
		{
			the_min_fit_inallG = min_fit;
			for (i = 0; i < INPUT_NEURAL_NUM; i++)
			{
				history_best_ch.nn_initial.input_w[i] = GA_best_ch.nn_initial.input_w[i];
				fprintf(f_history_ch, "%f\n", history_best_ch.nn_initial.input_w[i]);
			}	
			for (j = 0; j < SINGLE_WEIGHT; j++)
			{
				history_best_ch.nn_initial.w[j] = GA_best_ch.nn_initial.w[j];
				fprintf(f_history_ch, "%f\n", history_best_ch.nn_initial.w[j]);
			}
			history_best_ch.A = GA_best_ch.A;
			history_best_ch.B = GA_best_ch.B;
			history_best_ch.C = GA_best_ch.C;
			history_best_ch.D = GA_best_ch.D;
			fprintf(f_history_ch, "A:%f\n B:%f\n C:%f\n D:%f\n", history_best_ch.A, history_best_ch.B, history_best_ch.C, history_best_ch.D);
			record_g = g;
		}

		if (min_fit < the_min_fit_inallG)
		{
			the_min_fit_inallG = min_fit;
			for (i = 0; i < INPUT_NEURAL_NUM; i++)
			{
				history_best_ch.nn_initial.input_w[i] = GA_best_ch.nn_initial.input_w[i];
				fprintf(f_history_ch, "%f\n", history_best_ch.nn_initial.input_w[i]);
			}
			for (j = 0; j < SINGLE_WEIGHT; j++)
			{
				history_best_ch.nn_initial.w[j] = GA_best_ch.nn_initial.w[j];
				fprintf(f_history_ch, "%f\n", history_best_ch.nn_initial.w[j]);
			}
			history_best_ch.A = GA_best_ch.A;
			history_best_ch.B = GA_best_ch.B;
			history_best_ch.C = GA_best_ch.C;
			history_best_ch.D = GA_best_ch.D;
			fprintf(f_history_ch, "A:%f\n B:%f\n C:%f\n D:%f\n", history_best_ch.A, history_best_ch.B, history_best_ch.C, history_best_ch.D);
			record_g = g;
		}
		printf("the history min ranking is:%f \n", the_min_fit_inallG);
		fprintf(f_results, "history best is in g=%d:%f\n", record_g + 1, the_min_fit_inallG);
		cudaFree(d_fitness);
		cudaFree(d_ch);
		cudaFree(d_in_fitness);
		cudaFree(d_in_position);

	}//GA代
	fclose(f_history_ch);
	//fclose(f_position);
	fclose(f_results);
	free(h_ch);
	free(h_fitness);
	free(h_in_fitness);
	free(h_in_position);
	return 0;
}