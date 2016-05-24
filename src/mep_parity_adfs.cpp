// ---------------------------------------------------------------------------
// Multi Expression Programming - basic source code for solving Even Parity problems with Automatically Defined Functions
// (c) Mihai Oltean  www.mepx.org ; mihai.oltean@gmail.com
// Last update on: 2016.05.22

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// ---------------------------------------------------------------------------

// More info at:
// www.mepx.org
// www.mep.cs.ubbcluj.ro
// www.github.com/mepx

// paper to read:
// Oltean Mihai, Improving Multi Expression Programming: an Ascending Trail from Sea - level Even-3-parity Problem to Alpine Even-18-Parity Problem, chapter 10, Evolvable Machines: Theory and Applications, Springer - Verlag, edited by Nadia Nedjah(et al.), pp. 229 - 255, 2004
// http://www.mep.cs.ubbcluj.ro/oltean_parity.pdf

// Compiled with Microsoft Visual C++ 2013, XCode 7 and Embarcadero C++Builder XE

// Please reports any sugestions and/or bugs to mihai.oltean@gmail.com

// Training data file must have the following format (see even_7_parity.txt):

// m n
// x11 x12 ... x1n t1
// x21 x22 ....x2n t2
// .............
// xm1 xm2 ... xmn tm

// where m is the number of training data
// and n is the number of variables.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define max_num_inputs_ADF 4

#define num_operators 4 // 4 + 2 Automatically Defined Functions

// and -1
// or -2
// nand -3
// nor -4
// adf_0 -5
// adf_1 -6
// adf_2 -7

char operators_string[7][10] = {"and", "or", "nand", "nor", "adf_0", "adf_1", "adf_2"};

// ---------------------------------------------------------------------------
struct t_code3 {
	int op; // either a variable or an operator
	int addr1, addr2, addr3, addr4; // pointers to arguments
};
// ---------------------------------------------------------------------------
struct t_mep_chromosome {
	t_code3 *prg; // a string of genes
	t_code3 *adf_0; // ADF with 2 parameters
	t_code3 *adf_1; // ADF with 3 parameters
	t_code3 *adf_2; // ADF with 4 parameters
	int fitness;
	int best_index; // the index of the best expression in chromosome
};
// ---------------------------------------------------------------------------
struct t_mep_parameters {
	int code_length; // number of instructions in a chromosome
	int adf_code_length;
	int num_generations;
	int pop_size; // population size
	double mutation_probability, crossover_probability;
	double variables_probability, operators_probability;
	int num_adf_parameters[3];// for 3 ADF
};
// ---------------------------------------------------------------------------
void allocate_chromosome(t_mep_chromosome &c, t_mep_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	c.adf_0 = new t_code3[params.adf_code_length];
	c.adf_1 = new t_code3[params.adf_code_length];
	c.adf_2 = new t_code3[params.adf_code_length];
}
// ---------------------------------------------------------------------------
void delete_chromosome(t_mep_chromosome &c)
{
	if (c.prg) {
		delete[]c.prg;
		c.prg = NULL;
	}
	if (c.adf_0) {
		delete[]c.adf_0;
		c.adf_0 = NULL;
	}
	if (c.adf_1) {
		delete[]c.adf_1;
		c.adf_1 = NULL;
	}
	if (c.adf_2) {
		delete[]c.adf_2;
		c.adf_2 = NULL;
	}
}

// ---------------------------------------------------------------------------
void allocate_training_data(int **&data, int *&target, int num_training_data, int num_variables)
{
	target = new int[num_training_data];
	data = new int*[num_training_data];
	for (int i = 0; i < num_training_data; i++)
		data[i] = new int[num_variables];
}

// ---------------------------------------------------------------------------
void allocate_partial_expression_values(int **&expression_value, int **&adf_value, int num_training_data, int code_length, int adf_code_length)
{
	expression_value = new int*[code_length];
	for (int i = 0; i < code_length; i++)
		expression_value[i] = new int[num_training_data];

	adf_value = new int*[adf_code_length];
	for (int i = 0; i < adf_code_length; i++)
		adf_value[i] = new int[num_training_data];
}

// ---------------------------------------------------------------------------
void delete_partial_expression_values(int **&expression_value, int **&adf_value, int code_length, int adf_code_length)
{
	if (expression_value) {
		for (int i = 0; i < code_length; i++)
			delete[]expression_value[i];
		delete[]expression_value;
	}

	if (adf_value) {
		for (int i = 0; i < adf_code_length; i++)
			delete[]adf_value[i];
		delete[]adf_value;
	}
}

// ---------------------------------------------------------------------------
bool read_training_data(const char *filename, int **&training_data, int *&target, int &num_training_data, int &num_variables)
{
	FILE* f = fopen(filename, "r");
	if (!f)
		return false;

	fscanf(f, "%d%d", &num_training_data, &num_variables);

	allocate_training_data(training_data, target, num_training_data, num_variables);

	for (int i = 0; i < num_training_data; i++) {
		for (int j = 0; j < num_variables; j++)
			fscanf(f, "%d", &training_data[i][j]);
		fscanf(f, "%d", &target[i]);
	}

	fclose(f);
	return true;
}

// ---------------------------------------------------------------------------
void delete_data(int **&data, int *&target, int num_training_data)
{
	if (data)
		for (int i = 0; i < num_training_data; i++)
			delete[]data[i];
	delete[]data;
	delete[]target;
}

// ---------------------------------------------------------------------------
void copy_individual(t_mep_chromosome& dest, const t_mep_chromosome& source, t_mep_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.adf_code_length; i++) {
		dest.adf_0[i] = source.adf_0[i];
		dest.adf_1[i] = source.adf_1[i];
		dest.adf_2[i] = source.adf_2[i];
	}
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
}

// ---------------------------------------------------------------------------
void generate_random_chromosome(t_mep_chromosome &a, t_mep_parameters &params, int num_variables) // randomly initializes the individuals
{
	// on the first position we can have only a variable

	a.prg[0].op = rand() % num_variables;

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.code_length; i++) {
		double p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.prg[i].op = -rand() % (num_operators + 2) - 1; // an operator
		else
			a.prg[i].op = rand() % num_variables; // a variable

		a.prg[i].addr1 = rand() % i;
		a.prg[i].addr2 = rand() % i;
		a.prg[i].addr3 = rand() % i;
		a.prg[i].addr4 = rand() % i;
	}

	// adf_0
	a.adf_0[0].op = rand() % params.num_adf_parameters[0];

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.adf_code_length; i++) {
		double p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.adf_0[i].op = -rand() % num_operators - 1; // an operator
		else
			a.adf_0[i].op = rand() % params.num_adf_parameters[0]; // a variable

		a.adf_0[i].addr1 = rand() % i;
		a.adf_0[i].addr2 = rand() % i;
	}

	// adf_1
	a.adf_1[0].op = rand() % params.num_adf_parameters[1];

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.adf_code_length; i++) {
		double p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.adf_1[i].op = -rand() % num_operators - 1; // an operator
		else
			a.adf_1[i].op = rand() % params.num_adf_parameters[1]; // a variable

		a.adf_1[i].addr1 = rand() % i;
		a.adf_1[i].addr2 = rand() % i;
	}

	// adf_2
	a.adf_2[0].op = rand() % params.num_adf_parameters[2];

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.adf_code_length; i++) {
		double p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a.adf_2[i].op = -rand() % num_operators - 1; // an operator
		else
			a.adf_2[i].op = rand() % params.num_adf_parameters[2]; // a variable

		a.adf_2[i].addr1 = rand() % i;
		a.adf_2[i].addr2 = rand() % i;
	}

}

// ---------------------------------------------------------------------------
int compute_ADF(t_code3* adf_code, int adf_code_length, int num_training_data, int **adf_inputs, int* target, int **adf_cache_matrix, int &best_index)
{
	int best_error = num_training_data + 1;
	int sum_errors;

	for (int i = 0; i < adf_code_length; i++) { // read the t_chromosome from top to down
		sum_errors = 0; // and compute the fitness of each expression by dynamic programming

		switch (adf_code[i].op) {
		case -1: // and
			for (int k = 0; k < num_training_data; k++) {
				adf_cache_matrix[i][k] = adf_cache_matrix[adf_code[i].addr1][k] && adf_cache_matrix[adf_code[i].addr2][k];
				sum_errors += abs(adf_cache_matrix[i][k] - target[k]);
			}
			break;
		case -2: // or
			for (int k = 0; k < num_training_data; k++) {
				adf_cache_matrix[i][k] = adf_cache_matrix[adf_code[i].addr1][k] || adf_cache_matrix[adf_code[i].addr2][k];
				sum_errors += abs(adf_cache_matrix[i][k] - target[k]);
			}
			break;
		case -3: // nand
			for (int k = 0; k < num_training_data; k++) {
				adf_cache_matrix[i][k] = !(adf_cache_matrix[adf_code[i].addr1][k] && adf_cache_matrix[adf_code[i].addr2][k]);
				sum_errors += abs(adf_cache_matrix[i][k] - target[k]);
			}
			break;
		case -4: // nor
			for (int k = 0; k < num_training_data; k++) {
				adf_cache_matrix[i][k] = !(adf_cache_matrix[adf_code[i].addr1][k] || adf_cache_matrix[adf_code[i].addr2][k]);
				sum_errors += abs(adf_cache_matrix[i][k] - target[k]);
			}
			break;

		default: // a variable (parameter of adf)
			for (int k = 0; k < num_training_data; k++) {
				adf_cache_matrix[i][k] = adf_inputs[adf_code[i].op][k];
				sum_errors += abs(adf_cache_matrix[i][k] - target[k]);
			}
			break;
		}
		if (best_error > sum_errors) {
			best_error = sum_errors;
			best_index = i;
		}
	}
	return best_error;
}

// ---------------------------------------------------------------------------
void fitness(t_mep_chromosome &c, int code_length, int adf_code_length, int num_variables, int num_training_data, int **training_data, int *target, int **cache_matrix, int **adf_cache_matrix)
{
	// evaluate Individual
	c.fitness = num_training_data + 1;
	c.best_index = -1;
	int sum_errors;
	int *adf_inputs[max_num_inputs_ADF];
	int best_adf_index;

	for (int i = 0; i < code_length; i++) { // read the t_chromosome from top to down
		sum_errors = 0; // and compute the fitness of each expression by dynamic programming

		switch (c.prg[i].op) {
		case -1: // and
			for (int k = 0; k < num_training_data; k++) {
				cache_matrix[i][k] = cache_matrix[c.prg[i].addr1][k] && cache_matrix[c.prg[i].addr2][k];
				sum_errors += abs(cache_matrix[i][k] - target[k]);
			}
			break;
		case -2: // or
			for (int k = 0; k < num_training_data; k++) {
				cache_matrix[i][k] = cache_matrix[c.prg[i].addr1][k] || cache_matrix[c.prg[i].addr2][k];
				sum_errors += abs(cache_matrix[i][k] - target[k]);
			}
			break;
		case -3: // nand
			for (int k = 0; k < num_training_data; k++) {
				cache_matrix[i][k] = !(cache_matrix[c.prg[i].addr1][k] && cache_matrix[c.prg[i].addr2][k]);
				sum_errors += abs(cache_matrix[i][k] - target[k]);
			}
			break;
		case -4: // nor
			for (int k = 0; k < num_training_data; k++) {
				cache_matrix[i][k] = !(cache_matrix[c.prg[i].addr1][k] || cache_matrix[c.prg[i].addr2][k]);
				sum_errors += abs(cache_matrix[i][k] - target[k]);
			}
			break;
		case -5: // adf_0
			adf_inputs[0] = cache_matrix[c.prg[i].addr1];
			adf_inputs[1] = cache_matrix[c.prg[i].addr2];
			sum_errors = compute_ADF(c.adf_0, adf_code_length, num_training_data, adf_inputs, target, adf_cache_matrix, best_adf_index);
			for (int kk = 0; kk < num_training_data; kk++)
				cache_matrix[i][kk] = adf_cache_matrix[best_adf_index][kk];
			break;
		case -6: // adf_1
			adf_inputs[0] = cache_matrix[c.prg[i].addr1];
			adf_inputs[1] = cache_matrix[c.prg[i].addr2];
			adf_inputs[2] = cache_matrix[c.prg[i].addr3];
			sum_errors = compute_ADF(c.adf_1, adf_code_length, num_training_data, adf_inputs, target, adf_cache_matrix, best_adf_index);
			for (int kk = 0; kk < num_training_data; kk++)
				cache_matrix[i][kk] = adf_cache_matrix[best_adf_index][kk];
			break;

		case -7: // adf_2
			adf_inputs[0] = cache_matrix[c.prg[i].addr1];
			adf_inputs[1] = cache_matrix[c.prg[i].addr2];
			adf_inputs[2] = cache_matrix[c.prg[i].addr3];
			adf_inputs[3] = cache_matrix[c.prg[i].addr4];
			sum_errors = compute_ADF(c.adf_2, adf_code_length, num_training_data, adf_inputs, target, adf_cache_matrix, best_adf_index);
			for (int kk = 0; kk < num_training_data; kk++)
				cache_matrix[i][kk] = adf_cache_matrix[best_adf_index][kk];
			break;

		default: // a variable
			for (int k = 0; k < num_training_data; k++) {
				cache_matrix[i][k] = training_data[k][c.prg[i].op];
				sum_errors += abs(cache_matrix[i][k] - target[k]);
			}
			break;
		}

		if (c.fitness > sum_errors) {
			c.fitness = sum_errors;
			c.best_index = i;
		}
	}
}

// ---------------------------------------------------------------------------
void mutation(t_mep_chromosome &a_chromosome, t_mep_parameters params, int num_variables) // mutate the individual
{
	// main program
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability)
		a_chromosome.prg[0].op = rand() % num_variables;

	// other genes
	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX; // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -rand() % (num_operators + 2) - 1;
			else
				a_chromosome.prg[i].op = rand() % num_variables;
		}

		p = rand() / (double)RAND_MAX; // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr1 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr2 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the 3rd address   (addr3)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr3 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the 4th address   (addr4)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr4 = rand() % i;
	}

	// adf_0
	p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability)
		a_chromosome.adf_0[0].op = rand() % params.num_adf_parameters[0];

	// other genes
	for (int i = 1; i < params.adf_code_length; i++) {
		p = rand() / (double)RAND_MAX; // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.adf_0[i].op = -rand() % num_operators - 1;
			else
				a_chromosome.adf_0[i].op = rand() % params.num_adf_parameters[0];
		}

		p = rand() / (double)RAND_MAX; // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.adf_0[i].addr1 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.adf_0[i].addr2 = rand() % i;
	}

	// adf_1
	p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability)
		a_chromosome.adf_1[0].op = rand() % params.num_adf_parameters[1];

	// other genes
	for (int i = 1; i < params.adf_code_length; i++) {
		p = rand() / (double)RAND_MAX; // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.adf_1[i].op = -rand() % num_operators - 1;
			else
				a_chromosome.adf_1[i].op = rand() % params.num_adf_parameters[1];
		}

		p = rand() / (double)RAND_MAX; // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.adf_1[i].addr1 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.adf_1[i].addr2 = rand() % i;
	}

	// adf_2
	p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability)
		a_chromosome.adf_2[0].op = rand() % params.num_adf_parameters[2];

	// other genes
	for (int i = 1; i < params.adf_code_length; i++) {
		p = rand() / (double)RAND_MAX; // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.adf_2[i].op = -rand() % num_operators - 1;
			else
				a_chromosome.adf_2[i].op = rand() % params.num_adf_parameters[2];
		}

		p = rand() / (double)RAND_MAX; // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.adf_2[i].addr1 = rand() % i;

		p = rand() / (double)RAND_MAX; // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.adf_2[i].addr2 = rand() % i;
	}
}

// ---------------------------------------------------------------------------
void one_cut_point_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, t_mep_parameters &params, t_mep_chromosome &offspring1, t_mep_chromosome &offspring2)
{
	// main program
	int cutting_pct = rand() % params.code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.prg[i] = parent1.prg[i];
		offspring2.prg[i] = parent2.prg[i];
	}
	for (int i = cutting_pct; i < params.code_length; i++) {
		offspring1.prg[i] = parent2.prg[i];
		offspring2.prg[i] = parent1.prg[i];
	}

	// adf_0
	cutting_pct = rand() % params.adf_code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.adf_0[i] = parent1.adf_0[i];
		offspring2.adf_0[i] = parent2.adf_0[i];
	}
	for (int i = cutting_pct; i < params.adf_code_length; i++) {
		offspring1.adf_0[i] = parent2.adf_0[i];
		offspring2.adf_0[i] = parent1.adf_0[i];
	}
	// adf_1
	cutting_pct = rand() % params.adf_code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.adf_1[i] = parent1.adf_1[i];
		offspring2.adf_1[i] = parent2.adf_1[i];
	}
	for (int i = cutting_pct; i < params.adf_code_length; i++) {
		offspring1.adf_1[i] = parent2.adf_1[i];
		offspring2.adf_1[i] = parent1.adf_1[i];
	}

	// adf_2
	cutting_pct = rand() % params.adf_code_length;
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.adf_2[i] = parent1.adf_2[i];
		offspring2.adf_2[i] = parent2.adf_2[i];
	}
	for (int i = cutting_pct; i < params.adf_code_length; i++) {
		offspring1.adf_2[i] = parent2.adf_2[i];
		offspring2.adf_2[i] = parent1.adf_2[i];
	}
}

// ---------------------------------------------------------------------------
void uniform_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, t_mep_parameters &params, t_mep_chromosome &offspring1, t_mep_chromosome &offspring2)
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

	// adf_0
	for (int i = 0; i < params.adf_code_length; i++)
		if (rand() % 2) {
			offspring1.adf_0[i] = parent1.adf_0[i];
			offspring2.adf_0[i] = parent2.adf_0[i];
		}
		else {
			offspring1.adf_0[i] = parent2.adf_0[i];
			offspring2.adf_0[i] = parent1.adf_0[i];
		}

	// adf_1
	for (int i = 0; i < params.adf_code_length; i++)
		if (rand() % 2) {
			offspring1.adf_1[i] = parent1.adf_1[i];
			offspring2.adf_1[i] = parent2.adf_1[i];
		}
		else {
			offspring1.adf_1[i] = parent2.adf_1[i];
			offspring2.adf_1[i] = parent1.adf_1[i];
		}

		// adf_2
	for (int i = 0; i < params.adf_code_length; i++)
			if (rand() % 2) {
				offspring1.adf_2[i] = parent1.adf_2[i];
				offspring2.adf_2[i] = parent2.adf_2[i];
			}
			else {
				offspring1.adf_2[i] = parent2.adf_2[i];
				offspring2.adf_2[i] = parent1.adf_2[i];
			}

}

// ---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{ // comparator for quick sort
	if (((t_mep_chromosome*)a)->fitness > ((t_mep_chromosome*)b)->fitness)
		return 1;
	else if (((t_mep_chromosome*)a)->fitness < ((t_mep_chromosome*)b)->fitness)
		return -1;
	else
		return 0;
}

// ---------------------------------------------------------------------------
void print_chromosome(t_mep_chromosome& a, t_mep_parameters &params, int num_variables)
{
	printf("The best solution is:\n");

	printf("\nADF 0:\n");
	for (int i = 0; i < params.adf_code_length; i++)
		if (a.adf_0[i].op < 0)
			printf("%d: %s %d %d\n", i, operators_string[abs(a.adf_0[i].op) - 1], a.adf_0[i].addr1, a.adf_0[i].addr2);
		else
			printf("%d: ADF inputs[%d]\n", i, a.adf_0[i].op);

	printf("\nADF 1:\n");
	for (int i = 0; i < params.adf_code_length; i++)
		if (a.adf_1[i].op < 0)
			printf("%d: %s %d %d\n", i, operators_string[abs(a.adf_1[i].op) - 1], a.adf_1[i].addr1, a.adf_1[i].addr2);
		else
			printf("%d: ADF inputs[%d]\n", i, a.adf_1[i].op);

	printf("\nADF 2:\n");
	for (int i = 0; i < params.adf_code_length; i++)
		if (a.adf_2[i].op < 0)
			printf("%d: %s %d %d\n", i, operators_string[abs(a.adf_2[i].op) - 1], a.adf_2[i].addr1, a.adf_2[i].addr2);
		else
			printf("%d: ADF inputs[%d]\n", i, a.adf_2[i].op);

	printf("\nmain program:\n");
		for (int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0)
			if (a.prg[i].op >= -4) // binary operators
					printf("%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2);
			else 
				if (a.prg[i].op == -5) // ADF0
				printf("%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2);
				else
					if (a.prg[i].op == -6)// ADF 1
						printf("%d: %s %d %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2, a.prg[i].addr3);
					else
						if (a.prg[i].op == -7) // ADF 2
							printf("%d: %s %d %d %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2, a.prg[i].addr3, a.prg[i].addr4);
		else
			printf("%d: inputs[%d]\n", i, a.prg[i].op);

	printf("Best instruction index = %d\n", a.best_index);
	printf("Fitness (best error) = %d\n", a.fitness);
}

// ---------------------------------------------------------------------------
int tournament_selection(t_mep_chromosome *pop, int pop_size, int tournament_size) // Size is the size of the tournament
{
	int r, p;
	p = rand() % pop_size;
	for (int i = 1; i < tournament_size; i++) {
		r = rand() % pop_size;
		p = pop[r].fitness < pop[p].fitness ? r : p;
	}
	return p;
}

// ---------------------------------------------------------------------------
void start_steady_state_mep(t_mep_parameters &params, int **training_data, int* target, int num_training_data, int num_variables) // Steady-State
{
	// a steady state approach:
	// we work with 1 population
	// newly created individuals will replace the worst existing ones (only if they are better).

	// allocate memory
	t_mep_chromosome *population;
	population = new t_mep_chromosome[params.pop_size];
	for (int i = 0; i < params.pop_size; i++)
		allocate_chromosome(population[i], params);

	t_mep_chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params);
	allocate_chromosome(offspring2, params);

	int ** eval_matrix;
	int ** adf_eval_matrix;
	allocate_partial_expression_values(eval_matrix, adf_eval_matrix, num_training_data, params.code_length, params.adf_code_length);

	// initialize
	for (int i = 0; i < params.pop_size; i++) {
		generate_random_chromosome(population[i], params, num_variables);
		fitness(population[i], params.code_length, params.adf_code_length, num_variables, num_training_data, training_data, target, eval_matrix, adf_eval_matrix);

	}
	// sort ascendingly by fitness
	qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);

	printf("generation %d, best fitness = %d\n", 0, population[0].fitness);

	for (int g = 1; g < params.num_generations; g++) { // for each generation
		for (int k = 0; k < params.pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(population, params.pop_size, 2);
			int r2 = tournament_selection(population, params.pop_size, 2);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < params.crossover_probability)
				one_cut_point_crossover(population[r1], population[r2], params, offspring1, offspring2);
			else { // no crossover so the offspring are a copy of the parents
				copy_individual(offspring1, population[r1], params);
				copy_individual(offspring2, population[r2], params);
			}
			// mutate the result and compute fitness
			mutation(offspring1, params, num_variables);

			fitness(offspring1, params.code_length, params.adf_code_length, num_variables, num_training_data, training_data, target, eval_matrix, adf_eval_matrix);

			// mutate the other offspring and compute fitness
			mutation(offspring2, params, num_variables);

			fitness(offspring2, params.code_length, params.adf_code_length, num_variables, num_training_data, training_data, target, eval_matrix, adf_eval_matrix);

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
		printf("generation %d, best fitness = %d\n", g, population[0].fitness);
		if (population[0].fitness == 0)
			break; // stop when error 0 is found
	}
	// print best t_chromosome
	print_chromosome(population[0], params, num_variables);

	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int i = 0; i < params.pop_size; i++)
		delete_chromosome(population[i]);
	delete[]population;

	delete_data(training_data, target, num_training_data);

	delete_partial_expression_values(eval_matrix, adf_eval_matrix, params.code_length, params.adf_code_length);
}

// --------------------------------------------------------------------
int main(void) {
	t_mep_parameters params;

	params.pop_size = 800; // the number of individuals in population  (must be an even number!)
	params.code_length = 300;
	params.adf_code_length = 50;
	params.num_generations = 1000; // the number of generations
	params.mutation_probability = 0.01; // mutation probability
	params.crossover_probability = 0.9; // crossover probability

	params.variables_probability = 0.5;
	params.operators_probability = 0.5; // sum must be 1
	
	params.num_adf_parameters[0] = 2; // 1st ADF has 2 parameters
	params.num_adf_parameters[1] = 3; // 2nd ADF has 3 parameters
	params.num_adf_parameters[2] = 4; // 3rd ADF has 4 parameters

	int num_training_data, num_variables;
	int** training_data, *target;

	if (!read_training_data("datasets\\even_6_parity.txt", training_data, target, num_training_data, num_variables)) {
		printf("Cannot find input file! Please specify the full path!");
		getchar();
		return 1;
	}

	srand(1);
	start_steady_state_mep(params, training_data, target, num_training_data, num_variables);

	printf("Press enter ...");
	getchar();

	return 0;
}
// --------------------------------------------------------------------
