//---------------------------------------------------------------------------
//   Multi Expression Programming - basic source code
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2016.02.03

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------

//   More info at:  www.mep.cs.ubbcluj.ro

//   Compiled with Microsoft Visual C++ 2013
//   Also compiled with XCode 5.

//   Please reports any sugestions and/or bugs to       mihai.oltean@gmail.com

//   Training data file must have the following format (see building1.txt and cancer1.txt):
//   building1 and cancer1 data were taken from PROBEN1

//   x11 x12 ... x1n f1
//   x21 x22 ....x2n f2
//   .............
//   xm1 xm2 ... xmn fm

//   where m is the number of training data
//   and n is the number of variables.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define num_operators 4

// +   -1
// -   -2
// *   -3
// /   -4

char operators_string[5] = "+-*/";

//---------------------------------------------------------------------------
struct code3{
	int op;				// either a variable, operator or constant; 
	// variables are indexed from 0: 0,1,2,...; 
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int adr1, adr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct chromosome{
	code3 *prg;        // the program - a string of genes
	double *constants; // an array of constants

	double fitness;        // the fitness (or the error)
	// for regression is computed as sum of abs differences between target and obtained
	// for classification is computed as the number of incorrectly classified data
	int best_index;        // the index of the best expression in chromosome
};
//---------------------------------------------------------------------------
struct parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int pop_size;                // population size
	double mutation_probability, crossover_probability;
	int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

	int problem_type; //0 - regression, 1 - classification
	double classification_threshold; // for classification problems only
};
//---------------------------------------------------------------------------
void allocate_chromosome(chromosome &c, parameters &params)
{
	c.prg = new code3[params.code_length];
	if (params.num_constants)
		c.constants = new double[params.num_constants];
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(chromosome &c)
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
bool read_training_data(const char *filename, char list_separator, double **&training_data, double *&target, int &num_training_data, int &num_variables)
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
void copy_individual(chromosome& dest, const chromosome& source, parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(chromosome &a, parameters &params, int num_variables) // randomly initializes the individuals
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
		double p = rand() / (double)RAND_MAX;

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
void compute_eval_matrix(chromosome &c, int code_length, int num_variables, int num_training_data, double **training_data, double *target, double **eval_matrix)
{
	// we keep intermediate values in a matrix because when an error occurs (like division by 0) we mutate that gene into a variables.
	// in such case it is faster to have all intermediate results until current gene, so that we don't have to recompute them again.

	bool is_error_case;  // division by zero, other errors


	for (int i = 0; i < code_length; i++)   // read the chromosome from top to down
	{
		is_error_case = false;
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
// evaluate Individual
void fitness_regression(chromosome &c, int code_length, int num_variables, int num_training_data, double **training_data, double *target, double **eval_matrix)
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
void fitness_classification(chromosome &c, int code_length, int num_variables, int num_training_data, double **training_data, double *target, double **eval_matrix)
{
	c.fitness = 1e+308;
	c.best_index = -1;

	compute_eval_matrix(c, code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	for (int i = 0; i < code_length; i++) {   // read the chromosome from top to down
		int count_incorrect_classified = 0;
		for (int k = 0; k < num_training_data; k++)
			if (eval_matrix[i][k] < 0) // the program tells me that this data is in class 0
				count_incorrect_classified += target[k];
			else // the program tells me that this data is in class 1
				count_incorrect_classified += fabs(1 - target[k]);// difference between obtained and expected

		if (c.fitness > count_incorrect_classified) {
			c.fitness = count_incorrect_classified;
			c.best_index = i;
		}
	}
}
//---------------------------------------------------------------------------
void mutation(chromosome &a_chromosome, parameters params, int num_variables) // mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		double p = rand() / (double)RAND_MAX * sum;

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
void one_cut_point_crossover(const chromosome &parent1, const chromosome &parent2, parameters &params, chromosome &offspring1, chromosome &offspring2)
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
void uniform_crossover(const chromosome &parent1, const chromosome &parent2, parameters &params, chromosome &offspring1, chromosome &offspring2)
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
	if (((chromosome *)a)->fitness > ((chromosome *)b)->fitness)
		return 1;
	else
		if (((chromosome *)a)->fitness < ((chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
void print_chromosome(chromosome& a, parameters &params, int num_variables)
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
int tournament_selection(chromosome *pop, int pop_size, int tournament_size)     // Size is the size of the tournament
{
	int r, p;
	p = rand() % pop_size;
	for (int i = 1; i < tournament_size; i++) {
		r = rand() % pop_size;
		p = pop[r].fitness < pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
void start_steady_state_mep(parameters &params, double **training_data, double* target, int num_training_data, int num_variables)       // Steady-State 
{
	// a steady state approach:
	// we work with 1 population
	// newly created individuals will replace the worst existing ones (only if they are better).

	// allocate memory
	chromosome *population;
	population = new chromosome[params.pop_size];
	for (int i = 0; i < params.pop_size; i++)
		allocate_chromosome(population[i], params);

	chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params);
	allocate_chromosome(offspring2, params);

	double ** eval_matrix;
	allocate_partial_expression_values(eval_matrix, num_training_data, params.code_length);

	// initialize
	for (int i = 0; i < params.pop_size; i++) {
		generate_random_chromosome(population[i], params, num_variables);
		if (params.problem_type == 0)
			fitness_regression(population[i], params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
		else
			fitness_classification(population[i], params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);

	}
	// sort ascendingly by fitness
	qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);

	printf("generation %d, best fitness = %lf\n", 0, population[0].fitness);

	for (int g = 1; g < params.num_generations; g++) {// for each generation
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
			if (params.problem_type == 0)
				fitness_regression(offspring1, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
			else
				fitness_classification(offspring1, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
			// mutate the other offspring and compute fitness
			mutation(offspring2, params, num_variables);
			if (params.problem_type == 0)
				fitness_regression(offspring2, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);
			else
				fitness_classification(offspring2, params.code_length, num_variables, num_training_data, training_data, target, eval_matrix);

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
		printf("generation %d, best fitness = %lf\n", g, population[0].fitness);
	}
	// print best chromosome
	print_chromosome(population[0], params, num_variables);

	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int i = 0; i < params.pop_size; i++)
		delete_chromosome(population[i]);
	delete[] population;

	delete_data(training_data, target, num_training_data);

	delete_partial_expression_values(eval_matrix, params.code_length);
}
//--------------------------------------------------------------------
int main(void)
{
	parameters params;

	params.pop_size = 100;						    // the number of individuals in population  (must be an even number!)
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

	params.problem_type = 0;             //0 - regression, 1 - classification; DONT FORGET TO SET IT
	params.classification_threshold = 0; // only for classification problems

	int num_training_data, num_variables;
	double** training_data, *target;

	if (!read_training_data("building1.txt", ' ', training_data, target, num_training_data, num_variables)) {
		printf("Cannot find file! Please specify the full path!");
		getchar();
		return 1;
	}

	srand(0);
	start_steady_state_mep(params, training_data, target, num_training_data, num_variables);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------