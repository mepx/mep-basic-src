//---------------------------------------------------------------------------
//	Multi Expression Programming - basic source code for solving binary classification problems
//	author: Mihai Oltean  
//	mihai.oltean@gmail.com
//	Last update on: 2024.3.3.0

//	License: MIT
//---------------------------------------------------------------------------

//   More info at:  
//     https://mepx.org
//     https://mepx.github.io
//     https://github.com/mepx

//   Please reports any sugestions and/or bugs to mihai.oltean@gmail.com

//   Training data file must have the following format (see building1.txt and cancer1.txt from the dataset folder):
//   Note that building1 and cancer1 data were taken from PROBEN1

//   x_11 x_12 ... x_1n f_1
//   x_21 x_22 ....x_2n f_2
//   .............
//   x_m1 x_m2 ... x_mn f_m

//   where m is the number of training data
//   and n is the number of variables.
//   x_ij are the inputs
//   f_i are the outputs

//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//---------------------------------------------------------------------------
#define NUM_Operators 5

#define ADD_OP -1 // +
#define DIF_OP -2 // -
#define MUL_OP -3 // *
#define DIV_OP -4 // /
#define SIN_OP -5 // /

char operators_as_string[NUM_Operators + 1] = "+-*/s";
//---------------------------------------------------------------------------
struct t_code3{
	int op;				// either a variable, an operator or a constant
	// variables are indexed from 0: 0, 1, 2, ...
	// constants are indexed from num_variables; the first constant has index num_variables
	// operators are stored as negative numbers -1, -2, -3...
	int addr1, addr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_mep_chromosome{
	t_code3 *code;        // the program - a string of genes
	double *constants;   // an array of constants

	double fitness;        // the fitness (or the error)
	// for binary classification is computed as the number of incorrectly classified data
	int best_index;        // the index of the best expression in chromosome
};
//---------------------------------------------------------------------------
struct t_mep_parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int pop_size;                // population size
	double mutation_probability, crossover_probability;
	int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

	double classification_threshold; // for classification problems only
};
//---------------------------------------------------------------------------
void allocate_chromosome(t_mep_chromosome &c, const t_mep_parameters &params)
{
	c.code = new t_code3[params.code_length];// the code
	if (params.num_constants)
		c.constants = new double[params.num_constants];// constants
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(t_mep_chromosome &c)
{
	if (c.code) {
		delete[] c.code;
		c.code = NULL;
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
{// allocate memory for the matrix storing the output of each expression for each training data
	// this is allocated once and then reused, in order to reduce the number of allocations/deletions
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
void delete_training_data(double **&data, double *&target, int num_training_data)
{
	if (data)
		for (int i = 0; i < num_training_data; i++)
			delete[] data[i];
	delete[] data;
	delete[] target;
}
//---------------------------------------------------------------------------
void copy_individual(t_mep_chromosome& dest, const t_mep_chromosome& source, const t_mep_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.code[i] = source.code[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(t_mep_chromosome &a, const t_mep_parameters &params, int num_variables) 
// randomly initializes the individuals
{
	// generate constants first
	for (int c = 0; c < params.num_constants; c++)
		a.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;

	// on the first position we can have only a variable or a constant
	double sum = params.variables_probability + params.constants_probability;
	double p = rand() / (double)RAND_MAX * sum;

	if (p <= params.variables_probability)
		a.code[0].op = rand() % num_variables;
	else
		a.code[0].op = num_variables + rand() % params.num_constants;

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.code[i].op = -rand() % NUM_Operators - 1;        // an operator
		else {
			if (p <= params.operators_probability + params.variables_probability)
				a.code[i].op = rand() % num_variables;     // a variable
			else
				a.code[i].op = num_variables + rand() % params.num_constants; // index of a constant
		}
		a.code[i].addr1 = rand() % i;
		a.code[i].addr2 = rand() % i;
	}
}
//---------------------------------------------------------------------------
void compute_eval_matrix(const t_mep_chromosome &c, int code_length, int num_variables,
	int num_training_data, 	const double **training_data, 
	double **eval_matrix)
{
	// we keep intermediate values in a matrix because when an error occurs (like division by 0) we mutate that gene into a variables.
	// in such case it is faster to have all intermediate results until current gene, so that we don't have to recompute them again.

	for (int i = 0; i < code_length; i++){   // read the chromosome from top to down
		bool is_error_case = false;// division by zero, other errors
		switch (c.code[i].op) {

		case  ADD_OP:  // +
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.code[i].addr1][k] + eval_matrix[c.code[i].addr2][k];
			break;
		case  DIF_OP:  // -
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.code[i].addr1][k] - eval_matrix[c.code[i].addr2][k];

			break;
		case  MUL_OP:  // *
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = eval_matrix[c.code[i].addr1][k] * eval_matrix[c.code[i].addr2][k];
			break;
		case  DIV_OP:  //  /
			for (int k = 0; k < num_training_data; k++)
				if (fabs(eval_matrix[c.code[i].addr2][k]) < 1e-6) // test if it is too small
					is_error_case = true;
			if (is_error_case) {                                           // an division by zero error occured !!!
				c.code[i].op = rand() % num_variables;   // the gene is mutated into a terminal
				for (int k = 0; k < num_training_data; k++)
					eval_matrix[i][k] = training_data[k][c.code[i].op];
			}
			else    // normal execution....
				for (int k = 0; k < num_training_data; k++)
					eval_matrix[i][k] = eval_matrix[c.code[i].addr1][k] / eval_matrix[c.code[i].addr2][k];
			break;
		case  SIN_OP:  // sin
			for (int k = 0; k < num_training_data; k++)
				eval_matrix[i][k] = sin(eval_matrix[c.code[i].addr1][k]);
			break;
		default:  // a variable
			for (int k = 0; k < num_training_data; k++)
				if (c.code[i].op < num_variables)
					eval_matrix[i][k] = training_data[k][c.code[i].op];
				else
					eval_matrix[i][k] = c.constants[c.code[i].op - num_variables];
			break;
		}
	}
}
//---------------------------------------------------------------------------
void fitness_binary_classification(t_mep_chromosome &c, const t_mep_parameters& params, int num_variables,
	int num_training_data, const double **training_data, const double *target, 
	double **eval_matrix)
{
	// use a threshold
	// if is less than threshold, it belongs to class 0
	// otherwise belongs to class 1

	c.fitness = 1e+308;
	c.best_index = -1;

	compute_eval_matrix(c, params.code_length, num_variables, num_training_data, training_data, eval_matrix);// compute the output of each expression

	for (int i = 0; i < params.code_length; i++) {   // read the chromosome from top to down
		int count_incorrect_classified = 0;
		for (int k = 0; k < num_training_data; k++)
			if (eval_matrix[i][k] < params.classification_threshold) // the program tells me that this data is in class 0
				count_incorrect_classified += (int)target[k];
			else // the program tells me that this data is in class 1
				count_incorrect_classified += abs(1 - (int)target[k]);// difference between obtained and expected

		if (c.fitness > count_incorrect_classified) {
			c.fitness = count_incorrect_classified;
			c.best_index = i;
		}
	}
}
//---------------------------------------------------------------------------
void mutation(t_mep_chromosome &a_chromosome, const t_mep_parameters &params, int num_variables) 
// mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		p = rand() / (double)RAND_MAX * sum;

		if (p <= params.variables_probability)
			a_chromosome.code[0].op = rand() % num_variables;
		else
			a_chromosome.code[0].op = num_variables + rand() % params.num_constants;
	}
	// other genes
	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.code[i].op = -rand() % NUM_Operators - 1;
			else {
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.code[i].op = rand() % num_variables;
				else
					a_chromosome.code[i].op = num_variables + rand() % params.num_constants; // index of a constant
			}
		}

		p = rand() / (double)RAND_MAX;      // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.code[i].addr1 = rand() % i;

		p = rand() / (double)RAND_MAX;      // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.code[i].addr2 = rand() % i;
	}
	// mutate the constants
	for (int c = 0; c < params.num_constants; c++) {
		p = rand() / (double)RAND_MAX;
		if (p < params.mutation_probability)
			a_chromosome.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;
	}

}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2)
{
	int cutting_pct = rand() % params.code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.code[i] = parent1.code[i];
		offspring2.code[i] = parent2.code[i];
	}
	for (int i = cutting_pct; i < params.code_length; i++) {
		offspring1.code[i] = parent2.code[i];
		offspring2.code[i] = parent1.code[i];
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
void uniform_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2)
{
	for (int i = 0; i < params.code_length; i++)
		if (rand() % 2) {
			offspring1.code[i] = parent1.code[i];
			offspring2.code[i] = parent2.code[i];
		}
		else {
			offspring1.code[i] = parent2.code[i];
			offspring2.code[i] = parent1.code[i];
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
	if (((t_mep_chromosome *)a)->fitness > ((t_mep_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_mep_chromosome *)a)->fitness < ((t_mep_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
void print_chromosome(const t_mep_chromosome& a, const t_mep_parameters &params, int num_variables)
{
	printf("The chromosome is:\n");
	printf("//------------------------------------------\n");
	for (int i = 0; i < params.num_constants; i++)
		printf("constants[%d] = %lf\n", i, a.constants[i]);

	for (int i = 0; i < params.code_length; i++) {
		if (a.code[i].op < 0) {
			if (a.code[i].op == SIN_OP)
				printf("%d: sin %d\n", i, a.code[i].addr1);
			else// binary operators
				printf("%d: %c %d %d\n", i, operators_as_string[abs(a.code[i].op) - 1], a.code[i].addr1, a.code[i].addr2);
		}
		else {
			if (a.code[i].op < num_variables)
				printf("%d: inputs[%d]\n", i, a.code[i].op);
			else
				printf("%d: constants[%d]\n", i, a.code[i].op - num_variables);
		}
	}
	printf("//------------------------------------------\n");
	printf("Best index (output provider) = %d\n", a.best_index);
	printf("Fitness = %lf\n", a.fitness);
}
//---------------------------------------------------------------------------
int tournament_selection(const t_mep_chromosome *population, int pop_size, int tournament_size)     
// Size parameter is the size of the tournament
{
	int p;
	p = rand() % pop_size;
	for (int i = 1; i < tournament_size; i++) {
		int r = rand() % pop_size;
		p = population[r].fitness < population[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
void start_steady_state_mep(t_mep_parameters &params, 
	const double **training_data, const double* target, int num_training_data, int num_variables) 
{
	// a steady state approach:
	// we work with 1 (one) population
	// newly created individuals will replace the worst existing ones (only if the offspring are better).

	// allocate memory
	t_mep_chromosome *population;
	population = new t_mep_chromosome[params.pop_size];
	for (int i = 0; i < params.pop_size; i++)
		allocate_chromosome(population[i], params);

	t_mep_chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params);
	allocate_chromosome(offspring2, params);

	double ** eval_matrix;
	allocate_partial_expression_values(eval_matrix, num_training_data, params.code_length);

	// initialize
	for (int i = 0; i < params.pop_size; i++) {
		generate_random_chromosome(population[i], params, num_variables);
		fitness_binary_classification(population[i], params, num_variables, num_training_data, training_data, target, eval_matrix);
	}
	// sort ascendingly by fitness
	qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);

	printf("generation %d, best fitness = %lf\n", 0, population[0].fitness);

	for (int generation = 1; generation < params.num_generations; generation++) {// for each generation
		for (int k = 0; k < params.pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(population, params.pop_size, 2);
			int r2 = tournament_selection(population, params.pop_size, 2);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < params.crossover_probability)
				one_cut_point_crossover(population[r1], population[r2], params, offspring1, offspring2);
			else {// no crossover so the offspring are a copy of the parents
				copy_individual(offspring1, population[r1], params);
				copy_individual(offspring2, population[r2], params);
			}
			// mutate the result and compute fitness
			mutation(offspring1, params, num_variables);
			fitness_binary_classification(offspring1, params, num_variables, num_training_data, training_data, target, eval_matrix);
			// mutate the other offspring and compute fitness
			mutation(offspring2, params, num_variables);
			fitness_binary_classification(offspring2, params, num_variables, num_training_data, training_data, target, eval_matrix);

			// replace the worst in the population
			if (offspring1.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring1, params);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
			if (offspring2.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring2, params);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
		}
		printf("generation %d, best fitness = %lf\n", generation, population[0].fitness);
	}
	// print best chromosome which is always on position 0 because of sorting
	print_chromosome(population[0], params, num_variables);

	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int i = 0; i < params.pop_size; i++)
		delete_chromosome(population[i]);
	delete[] population;

	delete_partial_expression_values(eval_matrix, params.code_length);
}
//--------------------------------------------------------------------
bool get_next_field(char *start_sir, char list_separator, char* dest, int & size, int &skip_size)
{
	skip_size = 0;
	while (start_sir[skip_size] && (start_sir[skip_size] != '\n') && (start_sir[skip_size] == list_separator))
		skip_size++;// skip separator at the beginning

	size = 0;
	while (start_sir[skip_size + size] && (start_sir[skip_size + size] != list_separator) && (start_sir[skip_size + size] != '\n')) // run until a find a separator or end of line or new line char
		size++;

	if (!size || !start_sir[skip_size + size])
		return false;
	strncpy(dest, start_sir + skip_size, size);
	dest[size] = '\0';
	return true;
}
// ---------------------------------------------------------------------------
bool read_data(const char *filename, char list_separator, double **&data, double *&target, int &num_data, int &num_variables)
{
	FILE* f = fopen(filename, "r");
	if (!f) {
		num_data = 0;
		num_variables = 0;
		return false;
	}

	char *buf = new char[10000];
	char * start_buf = buf;
	// count the number of training data and the number of variables
	num_data = 0;
	while (fgets(buf, 10000, f)) {
		if (strlen(buf) > 1)
			num_data++;
		if (num_data == 1) {
			num_variables = 0;

			char tmp_str[10000];
			int size;
			int skip_size;
			bool result = get_next_field(buf, list_separator, tmp_str, size, skip_size);
			while (result) {
				buf = buf + size + 1 + skip_size;
				result = get_next_field(buf, list_separator, tmp_str, size, skip_size);
				num_variables++;
			}
		}
		buf = start_buf;
	}
	delete[] start_buf;
	num_variables--;
	rewind(f);

	allocate_training_data(data, target, num_data, num_variables);

	for (int i = 0; i < num_data; i++) {
		for (int j = 0; j < num_variables; j++)
			fscanf(f, "%lf", &data[i][j]);
		fscanf(f, "%lf", &target[i]);
	}
	fclose(f);

	return true;
}
//---------------------------------------------------------------------------
int main(void)
{
	t_mep_parameters params;

	params.pop_size = 100;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 200;					// the number of generations
	params.mutation_probability = 0.1;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 5; // use 5 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	params.classification_threshold = 0; // only for classification problems

	int num_training_data, num_variables;
	double** training_data, *target;

	if (!read_data("cancer1.txt", ' ', training_data, target, num_training_data, num_variables)) {
	
	//if (!read_data("datasets/building1.txt", ' ', training_data, target, num_training_data, num_variables)) {
		printf("Cannot find file! Please specify the full path!");
		getchar();
		return 1;
	}

	srand(0);// random seed

	clock_t start_time = clock();
	// run MEP
	start_steady_state_mep(params, (const double**)training_data, (const double*)target, num_training_data, num_variables);

	clock_t end_time = clock();

	double running_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

	printf("Running time = %lf seconds\n", running_time);

	delete_training_data(training_data, target, num_training_data);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------