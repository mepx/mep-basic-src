//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads
//   Author: Mihai Oltean  (mihai.oltean@gmail.com)
//   Version: 2021.11.25.0

//   License: MIT
//---------------------------------------------------------------------------
//   Each subpopulation is evolved in a separate thread
//   The subpopulations have a circular structure
//   At the end of each generation we move, from each subpopulation, some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations and threads)
//---------------------------------------------------------------------------
//   Compiled with Microsoft Visual C++ 2019
//   Also compiled with XCode 13.
//   Requires C++11 or newer (for thread support)

//---------------------------------------------------------------------------
//   How to use it: 
//   	just create a console project and copy-paste the content this file in the main file of the project
//---------------------------------------------------------------------------
//   More info at:  
//     https://mepx.org
//     https://mepx.github.io
//     https://github.com/mepx

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com
//---------------------------------------------------------------------------
//   Training data file must have the following format (see building1.txt and cancer1.txt):
//   building1 and cancer1 data are taken from PROBEN1

//   x11 x12 ... x1n f1
//   x21 x22 ....x2n f2
//   .............
//   xm1 xm2 ... xmn fm

//   where m is the number of training data
//   and n is the number of variables.

//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <thread>
#include <mutex>
//--------------------------------------------------------------------

#define num_operators 4

// +   -1
// -   -2
// *   -3
// /   -4

char operators_string[5] = "+-*/";

#define PROBLEM_REGRESSION 0
#define PROBLEM_BINARY_CLASSIFICATION 1

//---------------------------------------------------------------------------
struct t_code3{
	int op;				// either a variable, operator or constant; 
	// variables are indexed from 0: 0,1,2,...; 
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int adr1, adr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_chromosome{
	t_code3 *prg;        // the program - a string of genes
	double *constants; // an array of constants

	double fitness;        // the fitness (or the error)
	// for regression is computed as sum of abs differences between target and obtained
	// for classification is computed as the number of incorrectly classified data
	int best_instruction_index;        // the index of the best expression in chromosome
};
//---------------------------------------------------------------------------
struct t_parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int num_sub_populations;       // number of subpopulations
	int sub_population_size;                // subpopulation size
	double mutation_probability, crossover_probability;
	int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

	int problem_type; //regression or binary classification
	double classification_threshold; // for binary classification problems only

	int num_threads; // num threads. 
	//for best performances the number of subpopulations should be multiple of num_threads.
	// num_thread should no exceed the number of processor cores.
};
//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, const t_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	if (params.num_constants)
		c.constants = new double[params.num_constants];
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(t_chromosome &c)
{
	if (c.prg) {
		delete[] c.prg;
		c.prg = NULL;
	}
	if (c.constants) {
		delete[] c.constants;
		c.constants = NULL;
	}
}
//---------------------------------------------------------------------------
void allocate_training_data(double **&data, double *&target, int num_training_data, int num_variables)
{
	target = new double[num_training_data];
	data = new double*[num_training_data];
	for (int i = 0; i < num_training_data; i++)
		data[i] = new double[num_variables];
}
//---------------------------------------------------------------------------
void allocate_partial_expression_values(double ***&expression_value, int num_training_data, int code_length, int num_threads)
{
	// partial values are stored in a matrix of size: code_length x num_training_data
	// for each thread we have a separate matrix
	expression_value = new double**[num_threads];
	for (int t = 0; t < num_threads; t++) {
		expression_value[t] = new double*[code_length];
		for (int i = 0; i < code_length; i++)
			expression_value[t][i] = new double[num_training_data];
	}
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(double ***&expression_value, int code_length, int num_threads)
{
	if (expression_value) {
		for (int t = 0; t < num_threads; t++) {
			for (int i = 0; i < code_length; i++)
				delete[] expression_value[t][i];
			delete[] expression_value[t];
		}
		delete[] expression_value;
	}
}
//---------------------------------------------------------------------------
bool get_next_field(char *start_sir, char list_separator, char* dest, int & size)
{
	size = 0;
	while (start_sir[size] && (start_sir[size] != list_separator) && (start_sir[size] != '\n'))
		size++;
	if (!size && !start_sir[size])
		return false;
	strncpy(dest, start_sir, size);
	dest[size] = '\0';
	return true;
}
// ---------------------------------------------------------------------------
bool read_training_data(const char *filename, char list_separator, 
		double **&training_data, double *&target, int &num_training_data, int &num_variables)
{
	FILE* f = fopen(filename, "r");
	if (!f)
		return false;

	char *buf = new char[10000];
	char * start_buf = buf;
	// count the number of training data and the number of variables
	num_training_data = 0;
	while (fgets(buf, 10000, f)) {
		if (strlen(buf) > 1)
			num_training_data++;
		if (num_training_data == 1) {
			num_variables = 0;

			char tmp_str[10000];
			int size;
			bool result = get_next_field(buf, list_separator, tmp_str, size);
			while (result) {
				buf = buf + size + 1;
				result = get_next_field(buf, list_separator, tmp_str, size);
				num_variables++;
			}
		}
		buf = start_buf;
	}
	delete[] start_buf;
	num_variables--;
	rewind(f);

	allocate_training_data(training_data, target, num_training_data, num_variables);

	for (int i = 0; i < num_training_data; i++) {
		for (int j = 0; j < num_variables; j++)
			fscanf(f, "%lf", &training_data[i][j]);
		fscanf(f, "%lf", &target[i]);
	}
	fclose(f);
	return true;
}
//---------------------------------------------------------------------------
void delete_data(double **&data, double *&target, int num_training_data)
{
	if (data)
		for (int i = 0; i < num_training_data; i++)
			delete[] data[i];
	delete[] data;
	delete[] target;
}
//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, const t_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.best_instruction_index = source.best_instruction_index;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(t_chromosome &a, const t_parameters &params, int num_variables) // randomly initializes the individuals
{
	// generate constants first
	for (int c = 0; c < params.num_constants; c++)
		a.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;

	// on the first position we can have only a variable or a constant
	double sum = params.variables_probability + params.constants_probability;
	double p = rand() / (double)RAND_MAX * sum;

	if (p <= params.variables_probability)
		a.prg[0].op = rand() % num_variables;
	else
		a.prg[0].op = num_variables + rand() % params.num_constants;

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.prg[i].op = -rand() % num_operators - 1;        // an operator
		else
			if (p <= params.operators_probability + params.variables_probability)
				a.prg[i].op = rand() % num_variables;     // a variable
			else
				a.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant

		a.prg[i].adr1 = rand() % i;
		a.prg[i].adr2 = rand() % i;
	}
}
//---------------------------------------------------------------------------
void compute_eval_matrix(t_chromosome &c, 
			int code_length, int num_variables, int num_training_data, const double **training_data, const double *target, 
			double **eval_matrix)
{
	// we keep intermediate values in a matrix because when an error occurs (like division by 0) we mutate that gene into a variables.
	// in such case it is faster to have all intermediate results until current gene, so that we don't have to recompute them again.


	for (int i = 0; i < code_length; i++){   // read the chromosome from top to down
	
		bool is_error_case = false;// division by zero, other errors
		switch (c.prg[i].op) {

		case  -1:  // +
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] + eval_matrix[c.prg[i].adr2][k];
			break;
		case  -2:  // -
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] - eval_matrix[c.prg[i].adr2][k];

			break;
		case  -3:  // *
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] * eval_matrix[c.prg[i].adr2][k];
			break;
		case  -4:  //  /
			for (int k = 0; k < num_training_data; k++)
				if (fabs(eval_matrix[c.prg[i].adr2][k]) < 1e-6) // a small constant
					is_error_case = true;
			if (is_error_case) {                                           // an division by zero error occured !!!
				c.prg[i].op = rand() % num_variables;   // the gene is mutated into a terminal
				for (int k = 0; k < num_training_data; k++)
					eval_matrix[i][k] = training_data[k][c.prg[i].op];
			}
			else    // normal execution....
				for (int k = 0; k < num_training_data; k++)
					eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] / eval_matrix[c.prg[i].adr2][k];
			break;
		default:  // a variable
			for (int k = 0; k < num_training_data; k++)
				if (c.prg[i].op < num_variables)
					eval_matrix[i][k] = training_data[k][c.prg[i].op];
				else
					eval_matrix[i][k] = c.constants[c.prg[i].op - num_variables];
			break;
		}
	}
}
//---------------------------------------------------------------------------
// evaluate the chromosome c
void fitness_regression(t_chromosome &c, int code_length, int num_variables, 
		int num_training_data, const double **training_data, const double *target, 
		double **eval_matrix)
{
	c.fitness = 1e+308;
	c.best_instruction_index = -1;

	compute_eval_matrix(c, code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	for (int i = 0; i < code_length; i++) {   // read the chromosome from top to down
		double sum_of_errors = 0;
		for (int k = 0; k < num_training_data; k++)
			sum_of_errors += fabs(eval_matrix[i][k] - target[k]);// difference between obtained and expected

		if (c.fitness > sum_of_errors) {
			c.fitness = sum_of_errors;
			c.best_instruction_index = i;
		}
	}
}
//---------------------------------------------------------------------------
void fitness_binary_classification(t_chromosome &c, const t_parameters& params, int num_variables, int num_training_data, 
		const double **training_data, const double *target, 
		double **eval_matrix)
{
	c.fitness = 1e+308;
	c.best_instruction_index = -1;

	compute_eval_matrix(c, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	for (int i = 0; i < params.code_length; i++) {   // read the chromosome from top to down
		int count_incorrect_classified = 0;
		for (int k = 0; k < num_training_data; k++)
			if (eval_matrix[i][k] < params.classification_threshold) // the program tells me that this data is in class 0
				count_incorrect_classified += (int)target[k];
			else // the program tells me that this data is in class 1
				count_incorrect_classified += (int)fabs(1 - target[k]);// difference between obtained and expected

		if (c.fitness > count_incorrect_classified) {
			c.fitness = count_incorrect_classified;
			c.best_instruction_index = i;
		}
	}
}
//---------------------------------------------------------------------------
void mutation(t_chromosome &a_chromosome, const t_parameters& params, int num_variables) // mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		p = rand() / (double)RAND_MAX * sum;

		if (p <= params.variables_probability)
			a_chromosome.prg[0].op = rand() % num_variables;
		else
			a_chromosome.prg[0].op = num_variables + rand() % params.num_constants;
	}
	// other genes
	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -rand() % num_operators - 1;
			else {
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.prg[i].op = rand() % num_variables;
				else
					a_chromosome.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant
			}
		}

		p = rand() / (double)RAND_MAX;      // mutate the first address  (adr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr1 = rand() % i;

		p = rand() / (double)RAND_MAX;      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr2 = rand() % i;
	}
	// mutate the constants
	for (int c = 0; c < params.num_constants; c++) {
		p = rand() / (double)RAND_MAX; 
		if (p < params.mutation_probability)
			a_chromosome.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;
	}

}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_chromosome &parent1, const t_chromosome &parent2, 
		const t_parameters &params, 
		t_chromosome &offspring1, t_chromosome &offspring2)
{
	int cutting_pct = rand() % params.code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.prg[i] = parent1.prg[i];
		offspring2.prg[i] = parent2.prg[i];
	}
	for (int i = cutting_pct; i < params.code_length; i++) {
		offspring1.prg[i] = parent2.prg[i];
		offspring2.prg[i] = parent1.prg[i];
	}
	// now the constants
	if (params.num_constants) {
		cutting_pct = rand() % params.num_constants;
		for (int i = 0; i < cutting_pct; i++) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		for (int i = cutting_pct; i < params.num_constants; i++) {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
	}
}
//---------------------------------------------------------------------------
void uniform_crossover(const t_chromosome &parent1, const t_chromosome &parent2, 
		const t_parameters &params, 
	t_chromosome &offspring1, t_chromosome &offspring2)
{
	for (int i = 0; i < params.code_length; i++)
		if (rand() % 2) {
			offspring1.prg[i] = parent1.prg[i];
			offspring2.prg[i] = parent2.prg[i];
		}
		else {
			offspring1.prg[i] = parent2.prg[i];
			offspring2.prg[i] = parent1.prg[i];
		}

	// constants
	for (int i = 0; i < params.num_constants; i++)
		if (rand() % 2) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		else {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
}
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{// comparator for quick sort
	if (((t_chromosome *)a)->fitness > ((t_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_chromosome *)a)->fitness < ((t_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
void print_chromosome(const t_chromosome& a, const t_parameters &params, int num_variables)
{
	printf("The chromosome is:\n");

	for (int i = 0; i < params.num_constants; i++)
		printf("constants[%d] = %lf\n", i, a.constants[i]);

	for (int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0)
			printf("%d: %c %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].adr1, a.prg[i].adr2);
		else
			if (a.prg[i].op < num_variables)
				printf("%d: inputs[%d]\n", i, a.prg[i].op);
			else
				printf("%d: constants[%d]\n", i, a.prg[i].op - num_variables);

	printf("\nBest instruction index = %d\n", a.best_instruction_index);
	printf("Fitness = %lf\n", a.fitness);
}
//---------------------------------------------------------------------------
int tournament_selection(const t_chromosome *a_sub_pop, int sub_pop_size, int tournament_size)     // Size is the size of the tournament
{
	int p;
	p = rand() % sub_pop_size;
	for (int i = 1; i < tournament_size; i++) {
		int r = rand() % sub_pop_size;
		p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
void evolve_one_subpopulation(int *current_subpop_index, std::mutex* mutex, 
			t_chromosome ** sub_populations, int generation_index, 
			const t_parameters &params,
			const double **training_data, const double* target, int num_training_data, int num_variables, 
			double ** eval_matrix)
{
	int pop_index = 0;
	while (*current_subpop_index < params.num_sub_populations) {// still more subpopulations to evolve?

		while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
		mutex->unlock();

		// pop_index is the index of the subpopulation evolved by the current thread
		if (pop_index < params.num_sub_populations) {
			t_chromosome *a_sub_population = sub_populations[pop_index];

			t_chromosome offspring1, offspring2;
			allocate_chromosome(offspring1, params);
			allocate_chromosome(offspring2, params);

			if (generation_index == 0) {
				for (int i = 0; i < params.sub_population_size; i++) {
					generate_random_chromosome(a_sub_population[i], params, num_variables);
					if (params.problem_type == PROBLEM_REGRESSION)
						fitness_regression(a_sub_population[i], params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
					else
						fitness_binary_classification(a_sub_population[i], params, num_variables, num_training_data, training_data, target, eval_matrix);

				}
				// sort ascendingly by fitness inside this population
				qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
			}
			else // next generations
				for (int k = 0; k < params.sub_population_size; k += 2) {
					// we increase by 2 because at each step we create 2 offspring

					// choose the parents using binary tournament
					int r1 = tournament_selection(a_sub_population, params.sub_population_size, 2);
					int r2 = tournament_selection(a_sub_population, params.sub_population_size, 2);
					// crossover
					double p_0_1 = rand() / double(RAND_MAX); // a random number between 0 and 1
					if (p_0_1 < params.crossover_probability)
						one_cut_point_crossover(a_sub_population[r1], a_sub_population[r2], params, offspring1, offspring2);
					else {// no crossover so the offspring are a copy of the parents
						copy_individual(offspring1, a_sub_population[r1], params);
						copy_individual(offspring2, a_sub_population[r2], params);
					}
					// mutate the result and compute fitness
					mutation(offspring1, params, num_variables);
					if (params.problem_type == PROBLEM_REGRESSION)
						fitness_regression(offspring1, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
					else
						fitness_binary_classification(offspring1, params, num_variables, num_training_data, training_data, target, eval_matrix);
					// mutate the other offspring too
					mutation(offspring2, params, num_variables);
					if (params.problem_type == PROBLEM_REGRESSION)
						fitness_regression(offspring2, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
					else
						fitness_binary_classification(offspring2, params, num_variables, num_training_data, training_data, target, eval_matrix);

					// replace the worst in the population
					if (offspring1.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
						copy_individual(a_sub_population[params.sub_population_size - 1], offspring1, params);
						qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
					}
					if (offspring2.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
						copy_individual(a_sub_population[params.sub_population_size - 1], offspring2, params);
						qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
					}
				}

			delete_chromosome(offspring1);
			delete_chromosome(offspring2);
		}
	}
}
//---------------------------------------------------------------------------
void start_steady_state_mep(const t_parameters &params, 
				const double **training_data, const double* target, 
			int num_training_data, int num_variables)
{
	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory for all sub populations
	t_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation 
	}

	// allocate memory for 
	double *** eval_matrix;
	allocate_partial_expression_values(eval_matrix, num_training_data, params.code_length, params.num_threads);

	// an array of threads. Each sub population is evolved by a thread
	std::thread *mep_threads = new std::thread[params.num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads
	
	// initial population (generation 0)
	int current_subpop_index = 0;
	for (int t = 0; t < params.num_threads; t++)
		mep_threads[t] = std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, 0, params, training_data, target, num_training_data, num_variables, eval_matrix[t]);

	for (int t = 0; t < params.num_threads; t++) 
		mep_threads[t].join(); // wait for all threads to execute	

	// find the best individual from the entire population
	int best_individual_index = 0; // the index of the subpopulation containing the best invidual
	for (int p = 1; p < params.num_sub_populations; p++)
		if (sub_populations[p][0].fitness < sub_populations[best_individual_index][0].fitness)
			best_individual_index = p;

	printf("Generation = %d, Best fitness = %lf\n", 0, sub_populations[best_individual_index][0].fitness);
	
	// evolve for a fixed number of generations
	for (int generation = 1; generation < params.num_generations; generation++) { // for each generation

		current_subpop_index = 0;
		for (int t = 0; t < params.num_threads; t++)
			mep_threads[t] = std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, generation, params, training_data, target, num_training_data, num_variables, eval_matrix[t]);

		for (int t = 0; t < params.num_threads; t++) 
			mep_threads[t].join();

		// find the best individual
		best_individual_index = 0; // the index of the subpopulation containing the best invidual
		for (int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_index][0].fitness)
				best_individual_index = p;
		printf("Generation = %d, Best fitness = %lf\n", generation, sub_populations[best_individual_index][0].fitness);

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (int p = 0; p < params.num_sub_populations; p++) {
			int  k = rand() % params.sub_population_size;// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}
	}

	delete[] mep_threads;

		// print best chromosome

	print_chromosome(sub_populations[best_individual_index][0], params, num_variables);
	// free memory

	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

	delete_partial_expression_values(eval_matrix, params.code_length, params.num_threads);
}
//--------------------------------------------------------------------
int main(void)
{
	t_parameters params;
	params.num_sub_populations = 4;
	params.sub_population_size = 100;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 100;					// the number of generations
	params.mutation_probability = 0.1;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 3; // use 3 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	params.problem_type = PROBLEM_REGRESSION;             //DON'T FORGET TO SET IT
	params.classification_threshold = 0; // only for binary classification problems

	params.num_threads = 4;

	int num_training_data, num_variables;
	double** training_data, *target;

	char file_name[1000];
	strcpy(file_name, "building1.txt");

	if (!read_training_data(file_name, ' ', training_data, target, num_training_data, num_variables)) {
		printf("Cannot find %s file! Please specify the full path!", file_name);
		getchar();
		return 1;
	}

	srand(0); 

	printf("evolving...\n");
	start_steady_state_mep(params, (const double**)training_data, (const double*)target, num_training_data, num_variables);

	delete_data(training_data, target, num_training_data);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------
