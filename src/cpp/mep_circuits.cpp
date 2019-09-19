//---------------------------------------------------------------------------
//   Multi Expression Programming source code - for designing digital circuits
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2016.02.03

//   License: MIT
//---------------------------------------------------------------------------

//   paper to read:
//   Oltean Mihai, Grosan C., Evolving Digital Circuits using Multi Expression Programming, 
//   NASA/DoD Conference on Evolvable Hardware, 24-26 June, Seattle, 
//   Edited by R. Zebulum (et. al), pages 87-90, IEEE Press, NJ, 2004

//   Compiled with Microsoft Visual C++ 2013

//   More info at:  
//     www.mepx.org
//     https://mepx.github.io
//     https://github.com/mepx

//   Please reports any sugestions and/or bugs to       mihai.oltean@gmail.com

//   Training data file must have the following format (see 3x3_multiplier.txt):

//   num_training num_vars num_outputs
//   x11 x12 ... x1n o11 o12 ... o1k
//   x21 x22 ....x2n o21 o22 ... o2k
//   .............
//   xm1 xm2 ... xmn om om2 ... omk

//   where 
//   m is the number of training data
//   n is the number of inputs
//   k  is the number of outputs of the circuit


#include <stdio.h>
#include <stdlib.h>


#define num_operators 3 // also called gates

// a and b   -1
// a and not b   -2
// a xor b   -3

char operators_string[3][10] = {"a&b", "a&!b", "a^b"};

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
	code3 *prg;          // the circuit - a string of genes

	int fitness;         // the fitness (or the error)
	                     // number of incorrectly computed outputs
	
	int *best_index;     // the index of the best output in chromosome
	int num_gates;       // the actual number of utilized gates
};
//---------------------------------------------------------------------------
struct s_parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int pop_size;                // population size
	double mutation_probability, crossover_probability;
	double variables_probability, operators_probability;
};
//---------------------------------------------------------------------------
void allocate_chromosome(chromosome &a_chromosome, int code_length, int num_outputs)
{
	a_chromosome.prg = new code3[code_length];
	a_chromosome.best_index = new int[num_outputs];
}
//---------------------------------------------------------------------------
void delete_chromosome(chromosome &a_chromosome)
{
	if (a_chromosome.prg) {
		delete[] a_chromosome.prg;
		a_chromosome.prg = NULL;

		delete[] a_chromosome.best_index;
		a_chromosome.best_index = NULL;
	}
}
//---------------------------------------------------------------------------
void allocate_training_data(int **&input, int **&target, int num_training_data, int num_variables, int num_outputs)
{
	target = new int*[num_training_data];
	input = new int*[num_training_data];
	for (int i = 0; i < num_training_data; i++) {
		input[i] = new int[num_variables];
		target[i] = new int[num_outputs];
	}
}
//---------------------------------------------------------------------------
void allocate_partial_expression_values(int **&expression_value, int num_training_data, int code_length)
{
	expression_value = new int*[code_length];
	for (int i = 0; i < code_length; i++)
		expression_value[i] = new int[num_training_data];
}
//---------------------------------------------------------------------------
void allocate_error_matrix(int **&error_matrix, int code_length, int num_outputs)
{
	error_matrix = new int*[code_length];
	for (int i = 0; i < code_length; i++)
		error_matrix[i] = new int[num_outputs];
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(int **&expression_value, int code_length)
{
	if (expression_value) {
		for (int i = 0; i < code_length; i++)
			delete[] expression_value[i];
		delete[] expression_value;
	}
}
//---------------------------------------------------------------------------
void delete_error_matrix(int **&error_matrix, int code_length)
{
	if (error_matrix) {
		for (int i = 0; i < code_length; i++)
			delete[] error_matrix[i];
		delete[] error_matrix;
	}
}
//---------------------------------------------------------------------------
bool read_training_data(const char *filename, char list_separator, int **&training_data, int **&target, int &num_training_data, int &num_variables, int &num_outputs)
{
	FILE* f = fopen(filename, "r");
	if (!f)
		return false;

	// on the first line of the file we have 3 values: num_training_data, num_variables, num_outputs
	fscanf(f, "%d%d%d", &num_training_data, &num_variables, &num_outputs);

	allocate_training_data(training_data, target, num_training_data, num_variables, num_outputs);

	for (int i = 0; i < num_training_data; i++) {
		for (int j = 0; j < num_variables; j++)
			fscanf(f, "%d", &training_data[i][j]);
		for (int j = 0; j < num_outputs; j++)
			fscanf(f, "%d", &target[i][j]);
	}
	fclose(f);
	return true;
}
//---------------------------------------------------------------------------
void delete_data(int **&input, int **&target, int num_training_data)
{
	if (input) {
		for (int i = 0; i < num_training_data; i++) {
			delete[] input[i];
			delete[] target[i];
		}
		delete[] input;
		delete[] target;
	}
}
//---------------------------------------------------------------------------
void copy_individual(chromosome& dest, const chromosome& source, int code_length, int num_outputs)
{
	for (int i = 0; i < code_length; i++)
		dest.prg[i] = source.prg[i];

	for (int i = 0; i < num_outputs; i++)
		dest.best_index[i] = source.best_index[i];

	dest.fitness = source.fitness;
	dest.num_gates = source.num_gates;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(chromosome &a_chromosome, s_parameters &params, int num_variables)
// randomly initializes an individual
{
	// on the first position we can have only a variable

	a_chromosome.prg[0].op = rand() % num_variables;

	// for all other genes we put either an operator or variable
	for (int i = 1; i < params.code_length; i++) {
		double p = rand() / (double)RAND_MAX;

		if (p <= params.operators_probability)
			a_chromosome.prg[i].op = -rand() % num_operators - 1;        // an operator
		else
			a_chromosome.prg[i].op = rand() % num_variables;     // a variable

		a_chromosome.prg[i].adr1 = rand() % i;
		a_chromosome.prg[i].adr2 = rand() % i;
	}
}
//---------------------------------------------------------------------------
void mark_all(chromosome& a_chromosome, bool *marked, int current_index)
{
	// recusively mark all genes of a program
	if (!marked[current_index]) {
		if (a_chromosome.prg[current_index].op < 0) { // if it is an operator... mark all its arguments
			marked[current_index] = true;
			mark_all(a_chromosome, marked, a_chromosome.prg[current_index].adr1);
			mark_all(a_chromosome, marked, a_chromosome.prg[current_index].adr2);
		}
	}
}
//---------------------------------------------------------------------------
// evaluate Individual
void fitness(chromosome &a_chromosome, int code_length, int num_variables, int num_training_data, int num_outputs, int **training_data, int **target, int **eval_matrix, int **error_matrix)
{
	a_chromosome.fitness = num_outputs * num_training_data + 1; // the worst error we could have + 1

	// we keep intermediate values in a matrix

	for (int i = 0; i < code_length; i++)   // read the chromosome from top to down
	{
		for (int o = 0; o < num_outputs; o++)
			error_matrix[i][o] = 0;

		switch (a_chromosome.prg[i].op) {

		case  -1:  // a and b
			for (int k = 0; k < num_training_data; k++) {
				eval_matrix[i][k] = eval_matrix[a_chromosome.prg[i].adr1][k] & eval_matrix[a_chromosome.prg[i].adr2][k]; // the error 
				for (int o = 0; o < num_outputs; o++)
					error_matrix[i][o] += abs(eval_matrix[i][k] - target[k][o]);// error matrix for each gene and outut
			}
			break;
		case  -2:  // a and not b
			for (int k = 0; k < num_training_data; k++) {
				eval_matrix[i][k] = eval_matrix[a_chromosome.prg[i].adr1][k] & !eval_matrix[a_chromosome.prg[i].adr2][k];
				for (int o = 0; o < num_outputs; o++)
					error_matrix[i][o] += abs(eval_matrix[i][k] - target[k][o]);
			}
			break;
		case  -3:  //  a xor b
			for (int k = 0; k < num_training_data; k++) {
				eval_matrix[i][k] = eval_matrix[a_chromosome.prg[i].adr1][k] ^ eval_matrix[a_chromosome.prg[i].adr2][k];
				for (int o = 0; o < num_outputs; o++)
					error_matrix[i][o] += abs(eval_matrix[i][k] - target[k][o]);
			}
			break;
		default:  // a variable
			for (int k = 0; k < num_training_data; k++) {
				eval_matrix[i][k] = training_data[k][a_chromosome.prg[i].op];
				for (int o = 0; o < num_outputs; o++)
					error_matrix[i][o] += abs(eval_matrix[i][k] - target[k][o]);
			}
			break;
		}
	}
	// now I have to assign genes to outputs
	a_chromosome.fitness = 0;
	bool *used_gene = new bool[code_length];
	bool *used_output = new bool[num_outputs];
	bool *marked = new bool[code_length];

	for (int i = 0; i < code_length; marked[i++] = false);
	for (int i = 0; i < code_length; used_gene[i++] = false);
	for (int i = 0; i < num_outputs; used_output[i++] = false);

	for (int i = 0; i < num_outputs; i++) { // must repeat num_outputs because we have to find which genes provides which output
		int min_error = num_outputs * num_training_data + 1; // max error possible + 1
		int selected_gene = -1;
		int selected_output = -1;
		for (int g = 0; g < code_length; g++) // 
			if (!used_gene[g])
				for (int o = 0; o < num_outputs; o++)// find the min error not assigned yet
					if ((min_error > error_matrix[g][o]) && !used_output[o]) {
						min_error = error_matrix[g][o];
						selected_gene = g;
						selected_output = o;
					}
		used_gene[selected_gene] = true;
		used_output[selected_output] = true;
		mark_all(a_chromosome, marked, selected_gene);
		a_chromosome.fitness += min_error;
		a_chromosome.best_index[selected_output] = selected_gene;
	}
	// compute the number of gates used by the best circuit in the current chromosome
	a_chromosome.num_gates = 0;
	for (int i = 0; i < code_length; i++)
		if (marked[i])
			a_chromosome.num_gates++;

	delete[] used_gene;
	delete[] used_output;
	delete[] marked;

}
//---------------------------------------------------------------------------
void mutation(chromosome &a_chromosome, s_parameters params, int num_variables) // mutate the individual
{
	// mutate each symbol with the given probability
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability) 
		a_chromosome.prg[0].op = rand() % num_variables;

	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;      // mutate the operator
		if (p < params.mutation_probability) {
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -rand() % num_operators - 1;
			else
				a_chromosome.prg[i].op = rand() % num_variables;
		}

		p = rand() / (double)RAND_MAX;      // mutate the first address  (adr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr1 = rand() % i;

		p = rand() / (double)RAND_MAX;      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr2 = rand() % i;
	}
}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const chromosome &parent1, const chromosome &parent2, s_parameters &params, chromosome &offspring1, chromosome &offspring2)
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
}
//---------------------------------------------------------------------------
void uniform_crossover(const chromosome &parent1, const chromosome &parent2, s_parameters &params, chromosome &offspring1, chromosome &offspring2)
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
}
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{// comparator for quick sort
	chromosome * c1 = (chromosome *)a;
	chromosome * c2 = (chromosome *)b;
	if (c1->fitness > c2->fitness)
		return 1;
	else
		if (c1->fitness < c2->fitness)
			return -1;
		else
			if (c1->num_gates > c2->num_gates)
				return 1;
			else
				if (c1->num_gates < c2->num_gates)
					return -1;
				else
					return
					0;
}
//---------------------------------------------------------------------------
void print_chromosome(chromosome& a_chromosome, int code_length, int num_variables, int num_outputs)
{
	printf("The chromosome is:\n");

	for (int i = 0; i < code_length; i++)
		if (a_chromosome.prg[i].op < 0)
			printf("%d: %s %d %d\n", i, operators_string[abs(a_chromosome.prg[i].op) - 1], a_chromosome.prg[i].adr1, a_chromosome.prg[i].adr2);
		else
			printf("%d: inputs[%d]\n", i, a_chromosome.prg[i].op);

	for (int i = 0; i < num_outputs; i++)
		printf("output #%d is provided by gene #%d\n", i, a_chromosome.best_index[i]);

	printf("Fitness (error) = %d\n", a_chromosome.fitness);
	printf("Num gates = %d\n", a_chromosome.num_gates);
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
void start_steady_state_mep(s_parameters &params, int **training_data, int** target, int num_training_data, int num_variables, int num_outputs)       // Steady-State 
{
	// allocate memory
	chromosome *population;
	population = new chromosome[params.pop_size];
	for (int i = 0; i < params.pop_size; i++)
		allocate_chromosome(population[i], params.code_length, num_outputs);

	chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params.code_length, num_outputs);
	allocate_chromosome(offspring2, params.code_length, num_outputs);

	int ** eval_matrix;
	int **error_matrix;
	allocate_partial_expression_values(eval_matrix, num_training_data, params.code_length);
	allocate_error_matrix(error_matrix, params.code_length, num_outputs);

	// initialize
	for (int i = 0; i < params.pop_size; i++) {
		generate_random_chromosome(population[i], params, num_variables);
		fitness(population[i], params.code_length, num_variables, num_training_data, num_outputs, training_data, target, eval_matrix, error_matrix);
	}
	// sort ascendingly by fitness
	qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);

	printf("generation = %d, best error = %d, num_gates = %d\n", 0, population[0].fitness, population[0].num_gates);

	for (int g = 1; g < params.num_generations; g++) {
		for (int k = 0; k < params.pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(population, params.pop_size, 2);
			int r2 = tournament_selection(population, params.pop_size, 2);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < params.crossover_probability)
				uniform_crossover(population[r1], population[r2], params, offspring1, offspring2);
			else {// no crossover so the offspring are a copy of the parents
				copy_individual(offspring1, population[r1], params.code_length, num_outputs);
				copy_individual(offspring2, population[r2], params.code_length, num_outputs);
			}
			// mutate the result and move the mutant in the new population
			mutation(offspring1, params, num_variables);
			fitness(offspring1, params.code_length, num_variables, num_training_data, num_outputs, training_data, target, eval_matrix, error_matrix);

			mutation(offspring2, params, num_variables);
			fitness(offspring2, params.code_length, num_variables, num_training_data, num_outputs, training_data, target, eval_matrix, error_matrix);

			// replace the worst in the population
			if (offspring1.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring1, params.code_length, num_outputs);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
			if (offspring2.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring2, params.code_length, num_outputs);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
		}
		printf("generation = %d, best error = %d, num_gates = %d\n", g, population[0].fitness, population[0].num_gates);
	}

	fitness(population[0], params.code_length, num_variables, num_training_data, num_outputs, training_data, target, eval_matrix, error_matrix);
	// print best chromosome
	if (population[0].fitness == 0)
		print_chromosome(population[0], params.code_length, num_variables, num_outputs);
	else
		printf("No solution found! Please modify the parameters.\n");

	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int i = 0; i < params.pop_size; i++)
		delete_chromosome(population[i]);
	delete[] population;

	delete_data(training_data, target, num_training_data);

	delete_partial_expression_values(eval_matrix, params.code_length);
	delete_error_matrix(error_matrix, params.code_length);
}
//--------------------------------------------------------------------
int main(void)
{
	s_parameters params;

	params.pop_size = 200;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 30;						// max number of gates in the circuit
	params.num_generations = 1000;					// the number of generations
	params.mutation_probability = 0.01;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.5;
	params.operators_probability = 0.5;

	int num_training_data, num_variables, num_outputs;
	int** training_data, **target;

	if (!read_training_data("2x2_multiplier.txt", ' ', training_data, target, num_training_data, num_variables, num_outputs)) {
		printf("Cannot find file 2x2_multiplier.txt! Please specify the full path!");
		getchar();
		return 1;
	}

	srand(1);
	start_steady_state_mep(params, training_data, target, num_training_data, num_variables, num_outputs);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------