//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations
//   Author: Mihai Oltean (mihai.oltean@gmail.com)
//   Version: 2021.11.27

//   License: MIT
//---------------------------------------------------------------------------
//   The subpopulations have a circular structure
//   From each subpopulation we move some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations)
//---------------------------------------------------------------------------

//   Compiled with Microsoft Visual C++ 2019
//   Also compiled with XCode 5.
//---------------------------------------------------------------------------

//   How to use it: 
//   	just create a console project and copy-paste the content this file in the main file of the project

//   More info at:  
//     https://mepx.org
//     https://mepx.github.io
//     https://github.com/mepx

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com
//---------------------------------------------------------------------------

//   Training data file must have the following format (see building1.txt and cancer1.txt from datasets folder):
//   building1 and cancer1 data were taken from PROBEN1

//   x_11 x_12 ... x_1n f_1
//   x_21 x_22 ....x_2n f_2
//   .............
//   x_m1 x_m2 ... x_mn f_m

//   where m is the number of training data
//   and n is the number of variables.

//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PROBLEM_REGRESSION 0
#define PROBLEM_BINARY_CLASSIFICATION 1

#define num_operators 4

#define ADD_OP -1 // +
#define DIF_OP -2 // -
#define MUL_OP -3 // *
#define DIV_OP -4 // /


char operators_string[5] = "+-*/";

//---------------------------------------------------------------------------
struct t_code3{
	int op;				// either a variable, an operator or a constant 
	// variables are indexed from 0: 0, 1, 2, ...;
	// constants are indexed starting with num_variables index
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
	int best_index;        // the index of the best expression in chromosome
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

	int problem_type; //0 - regression, 1 - classification
	double classification_threshold; // for classification problems only
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
void allocate_partial_expression_values(double **&expression_value, int num_training_data, int code_length)
{
	expression_value = new double*[code_length];
	for (int i = 0; i < code_length; i++)
		expression_value[i] = new double[num_training_data];
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(double **&expression_value, int code_length)
{
	if (expression_value) {
		for (int i = 0; i < code_length; i++)
			delete[] expression_value[i];
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
void delete_training_data(double **&data, double *&target, int num_training_data)
{
	if (data) {
		for (int i = 0; i < num_training_data; i++)
			delete[] data[i];
		delete[] data;
		data = NULL;
	}
	if (target) {
		delete[] target;
		target = NULL;
	}
}
//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, const t_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
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
void compute_eval_matrix(t_chromosome &c, int code_length, int num_variables, int num_training_data, double **training_data, double *target, double **eval_matrix)
{
	// we keep intermediate values in a matrix because when an error occurs (like division by 0) we mutate that gene into a variables.
	// in such case it is faster to have all intermediate results until current gene, so that we don't have to recompute them again.

	for (int i = 0; i < code_length; i++){   // read the chromosome from top to down
		bool is_error_case = false;// division by zero, other errors
		switch (c.prg[i].op) {

		case  ADD_OP:  // +
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] + eval_matrix[c.prg[i].adr2][k];
			break;
		case  DIF_OP:  // -
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] - eval_matrix[c.prg[i].adr2][k];

			break;
		case  MUL_OP:  // *
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.prg[i].adr1][k] * eval_matrix[c.prg[i].adr2][k];
			break;
		case  DIV_OP:  //  /
			for (int k = 0; k < num_training_data; k++)
				if (fabs(eval_matrix[c.prg[i].adr2][k]) < 1e-6) // avoid division by near zero values
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
// evaluate Individual
void fitness_regression(t_chromosome &c, int code_length, int num_variables, int num_training_data, 
		double **training_data, double *target, double **eval_matrix)
{
	c.fitness = 1e+308;
	c.best_index = -1;

	compute_eval_matrix(c, code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	for (int i = 0; i < code_length; i++) {   // read the chromosome from top to down
		double sum_of_errors = 0;
		for (int k = 0; k < num_training_data; k++)
			sum_of_errors += fabs(eval_matrix[i][k] - target[k]);// difference between obtained and expected

		if (c.fitness > sum_of_errors) {
			c.fitness = sum_of_errors;
			c.best_index = i;
		}
	}
}
//---------------------------------------------------------------------------
void fitness_classification(t_chromosome &c, int code_length, int num_variables, int num_training_data, 
		double **training_data, double *target, double **eval_matrix)
{
	c.fitness = 1e+308;
	c.best_index = -1;

	compute_eval_matrix(c, code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	for (int i = 0; i < code_length; i++) {   // read the chromosome from top to down
		int count_incorrect_classified = 0;
		for (int k = 0; k < num_training_data; k++)
			if (eval_matrix[i][k] < 0) // the program tells me that this data is in class 0
				count_incorrect_classified += (int)target[k];
			else // the program tells me that this data is in class 1
				count_incorrect_classified += (int)fabs(1 - target[k]);// difference between obtained and expected

		if (c.fitness > count_incorrect_classified) {
			c.fitness = count_incorrect_classified;
			c.best_index = i;
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
			else
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.prg[i].op = rand() % num_variables;
				else
					a_chromosome.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant
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
void one_cut_point_crossover(
	const t_chromosome &parent1, const t_chromosome &parent2, 
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
void uniform_crossover(
	const t_chromosome &parent1, const t_chromosome &parent2, 
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

	printf("best index = %d\n", a.best_index);
	printf("Fitness = %lf\n", a.fitness);
}
//---------------------------------------------------------------------------
int tournament_selection(const t_chromosome *a_sub_pop, int sub_pop_size, int tournament_size)     // Size is the size of the tournament
{
	int p = rand() % sub_pop_size;
	for (int i = 1; i < tournament_size; i++) {
		int r = rand() % sub_pop_size;
		p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
void evolve_one_subpopulation(t_chromosome * a_sub_population, const t_parameters &params, 
			double **training_data, double* target, int num_training_data, int num_variables, 
			t_chromosome &offspring1, t_chromosome &offspring2, double ** eval_matrix)
{
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
			fitness_classification(offspring1, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
		// mutate the other offspring too
		mutation(offspring2, params, num_variables);
		if (params.problem_type == PROBLEM_REGRESSION)
			fitness_regression(offspring2, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
		else
			fitness_classification(offspring2, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);

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

}
//---------------------------------------------------------------------------
void start_steady_state_mep(const t_parameters &params, double **training_data, double* target, 
		int num_training_data, int num_variables)       // Steady-State 
{
	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory
	t_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation 
	}

	t_chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params);
	allocate_chromosome(offspring2, params);

	double ** eval_matrix;
	allocate_partial_expression_values(eval_matrix, num_training_data, params.code_length);

	// initialize
	for (int p = 0; p < params.num_sub_populations; p++)
		for (int i = 0; i < params.sub_population_size; i++) {
			generate_random_chromosome(sub_populations[p][i], params, num_variables);
			if (params.problem_type == PROBLEM_REGRESSION)
				fitness_regression(sub_populations[p][i], params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
			else
				fitness_classification(sub_populations[p][i], params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);

		}
	// sort ascendingly by fitness inside each population
	for (int p = 0; p < params.num_sub_populations; p++)
		qsort((void *)sub_populations[p], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);

	// find the best individual from the entire population
	int best_individual_index = 0; // the index of the subpopulation containing the best invidual
	for (int p = 1; p < params.num_sub_populations; p++)
		if (sub_populations[p][0].fitness < sub_populations[best_individual_index][0].fitness)
			best_individual_index = p;

	printf("generation %d, best fitness = %lf\n", 0, sub_populations[best_individual_index][0].fitness);

	// evolve for a fixed number of generations
	for (int g = 1; g < params.num_generations; g++) { // for each generation
		for (int p = 0; p < params.num_sub_populations; p++)  // for each subpopulation
			evolve_one_subpopulation(sub_populations[p], params, training_data, target, num_training_data, num_variables, offspring1, offspring2, eval_matrix);

		// find the best individual
		best_individual_index = 0; // the index of the subpopulation containing the best invidual
		for (int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_index][0].fitness)
				best_individual_index = p;
		printf("generation %d, best fitness = %lf\n", g, sub_populations[best_individual_index][0].fitness);

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		if (params.num_sub_populations > 1) // only if we have more than one sub_population
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
		// print best chromosome

	print_chromosome(sub_populations[best_individual_index][0], params, num_variables);
	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

	delete_partial_expression_values(eval_matrix, params.code_length);
}
//--------------------------------------------------------------------
int main(void)
{
	t_parameters params;
	params.num_sub_populations = 2;
	params.sub_population_size = 50;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 30;
	params.num_generations = 10;					// the number of generations
	params.mutation_probability = 0.1;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 3; // use 3 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	params.problem_type = PROBLEM_REGRESSION;             //0 - regression, 1 - classification; DONT FORGET TO SET IT
	params.classification_threshold = 0; // only for classification problems

	int num_training_data, num_variables;
	double** training_data, *target;

	if (!read_training_data("building1.txt", ' ', training_data, target, num_training_data, num_variables)) {
		printf("Cannot find building1.txt file! Please specify the full path!");
		getchar();
		return 1;
	}

	srand(0);
	printf("evolving...\n");
	start_steady_state_mep(params, training_data, target, num_training_data, num_variables);

	delete_training_data(training_data, target, num_training_data);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------
