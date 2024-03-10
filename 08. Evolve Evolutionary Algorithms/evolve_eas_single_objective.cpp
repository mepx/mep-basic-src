//---------------------------------------------------------------------------
//   Multi Expression Programming for evolving Evolutionary Algorithms
//   Author: Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2024.08.10.0
//   New versions of this program will be available at: https://github.com/mepx/mep-basic-src

//   Just create a console application and set this file as the main file of the project

//   MIT License

//   Paper to read:

//   Oltean Mihai, Grosan C., 
//   Evolving Evolutionary Algorithms using Multi Expression Programming, 
//   The 7th European Conference on Artificial Life, 
//   Dortmund, Edited by W.Banzhaf(et al), LNAI 2801, pp. 651 - 658, Springer - Verlag, Berlin, 2003.
//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
//---------------------------------------------------------------------------
#define NUM_MICRO_EA_OPERATORS 4

#define MICRO_EA_RANDOM_INIT 0
#define MICRO_EA_BINARY_SELECTION 1
#define MICRO_EA_CROSSOVER 2
#define MICRO_EA_MUTATION 3
//---------------------------------------------------------------------------
struct t_code3{// three address code
	int op;		// operators are the MICRO EA OPERATORS
	int adr1, adr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_meta_gp_chromosome{
	t_code3 *prg;        // the program - a string of genes

	double fitness;        // the fitness (or the error)
	int best_index;        // the index of the best expression in chromosome
};
//---------------------------------------------------------------------------
struct t_meta_gp_parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int pop_size;                // population size
	double mutation_probability, crossover_probability;
};
//---------------------------------------------------------------------------
struct t_micro_ea_parameters{
	double mutation_probability;
	int num_runs;
	int num_bits_per_dimension;
};
//---------------------------------------------------------------------------
void allocate_meta_chromosome(t_meta_gp_chromosome &c, t_meta_gp_parameters &params)
{
	c.prg = new t_code3[params.code_length];
}
//---------------------------------------------------------------------------
void delete_meta_chromosome(t_meta_gp_chromosome &c)
{
	if (c.prg) {
		delete[] c.prg;
		c.prg = NULL;
	}
}
//---------------------------------------------------------------------------
void copy_individual(t_meta_gp_chromosome& dest, const t_meta_gp_chromosome& source, t_meta_gp_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
}
//---------------------------------------------------------------------------
void generate_random_meta_chromosome(t_meta_gp_chromosome &a, t_meta_gp_parameters &params) // randomly initializes the individuals
{
	a.prg[0].op = MICRO_EA_RANDOM_INIT; // init only

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.code_length; i++) {
		a.prg[i].op = rand() % NUM_MICRO_EA_OPERATORS;

		a.prg[i].adr1 = rand() % i;
		a.prg[i].adr2 = rand() % i;
	}
}
//---------------------------------------------------------------------------
void mutate_meta_chromosome(t_meta_gp_chromosome &a_chromosome, t_meta_gp_parameters& params) // mutate the individual
{
	// mutate each symbol with the given probability
	// no mutation for the first gene
	// other genes

	for (int i = 1; i < params.code_length; i++) {
		double p = rand() / (double)RAND_MAX;      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = rand() / (double)RAND_MAX;

			a_chromosome.prg[i].op = rand() % NUM_MICRO_EA_OPERATORS;
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
void one_cut_point_crossover(const t_meta_gp_chromosome &parent1, const t_meta_gp_chromosome &parent2, t_meta_gp_parameters &params, t_meta_gp_chromosome &offspring1, t_meta_gp_chromosome &offspring2)
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
void uniform_crossover(const t_meta_gp_chromosome &parent1, const t_meta_gp_chromosome &parent2, t_meta_gp_parameters &params, t_meta_gp_chromosome &offspring1, t_meta_gp_chromosome &offspring2)
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
	if (((t_meta_gp_chromosome *)a)->fitness > ((t_meta_gp_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_meta_gp_chromosome *)a)->fitness < ((t_meta_gp_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
void print_meta_chromosome(const t_meta_gp_chromosome& an_individual, int code_length)
{
	printf("The chromosome is:\n");

	for (int i = 0; i < code_length; i++)
		switch (an_individual.prg[i].op) {
		case MICRO_EA_RANDOM_INIT:
			printf("%d: MICRO_EA_RANDOM_INIT\n", i);
			break;
		case MICRO_EA_BINARY_SELECTION:
			printf("%d: MICRO_EA_BINARY_SELECT(%d, %d)\n", i, an_individual.prg[i].adr1, an_individual.prg[i].adr2);
			break;
		case MICRO_EA_CROSSOVER:
			printf("%d: MICRO_EA_CROSSOVER(%d)\n", i, an_individual.prg[i].adr1);
			break;
		case MICRO_EA_MUTATION:
			printf("%d: MICRO_EA_MUTATION(%d, %d)\n", i, an_individual.prg[i].adr1, an_individual.prg[i].adr2);
			break;
	}
		
	printf("best index = %d\n", an_individual.best_index);
	printf("Fitness = %lf\n", an_individual.fitness);
}
//---------------------------------------------------------------------------
int tournament_selection(t_meta_gp_chromosome *pop, int pop_size, int tournament_size)     // Size is the size of the tournament
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
// function to be optimized
//---------------------------------------------------------------------------
// De Jong's function
// domeniul -5.12<=x(i)<=5.12
// minim f(x) = 0
// x = (0,0,..0)
double f(double *x, int num_dimensions)
{
	double result = 0;
	for (int i = 0; i < num_dimensions; i++)
		result = result + i * x[i] * x[i];
	return result;
}
//---------------------------------------------------------------------------
double binary_to_real(char* b_string, int num_bits_per_dimension, double min_x, double max_x)
{
	// transform a binary string of num_bits_per_dimension size into a real number in [min_x ... max_x] interval
	double x_real = 0;
	for (int j = 0; j < num_bits_per_dimension; j++)
		x_real = x_real * 2 + (int)b_string[j];    // now I have them in interval [0 ... 2 ^ num_bits_per_dimension - 1]
	x_real /= ((1 << num_bits_per_dimension) - 1); // now I have them in [0 ... 1] interval
	x_real *= (max_x - min_x);  // now I have them in [0 ... max_x - min_x] interval
	x_real += min_x;            // now I have them in [min_x ... max_x] interval

	return x_real;
}
//--------------------------------------------------------------------
void compute_fitness(t_meta_gp_chromosome &an_individual, int code_length, t_micro_ea_parameters &micro_params, int num_dimensions, double min_x, double max_x)
{
	double *x = new double[num_dimensions]; // buffer for storing real values
	double *micro_fitness = new double[code_length]; // fitness for each micro EA chromosome
	double *average_micro_fitness = new double[code_length]; // fitness for each micro EA chromosome

	// allocate some memory
	char **micro_values;
	micro_values = new char*[code_length];
	for (int i = 0; i < code_length; i++)
		micro_values[i] = new char[num_dimensions * micro_params.num_bits_per_dimension];

	for (int i = 0; i < code_length; i++)
		average_micro_fitness[i] = 0;
	
	// evaluate code
	for (int r = 0; r < micro_params.num_runs; r++) {// micro ea is run on multi runs

		for (int i = 0; i < code_length; i++) {
			switch (an_individual.prg[i].op) {
			case  MICRO_EA_RANDOM_INIT:       // Initialization
				for (int j = 0; j < num_dimensions; j++)
					for (int k = 0; k < micro_params.num_bits_per_dimension; k++)
						micro_values[i][j * micro_params.num_bits_per_dimension + k] = rand() % 2; // random values

				// compute fitness of that micro chromosome
				// transform to base 10
				for (int j = 0; j < num_dimensions; j++)
					x[j] = binary_to_real(micro_values[i], micro_params.num_bits_per_dimension, min_x, max_x);

				micro_fitness[i] = f(x, num_dimensions);// apply f - compute fitness of micro
				break;

			case MICRO_EA_BINARY_SELECTION:  // Selection (binary tournament)
				if (micro_fitness[an_individual.prg[i].adr1] < micro_fitness[an_individual.prg[i].adr2]) {
					memcpy(micro_values[i], micro_values[an_individual.prg[i].adr1], micro_params.num_bits_per_dimension * num_dimensions);
					micro_fitness[i] = micro_fitness[an_individual.prg[i].adr1];
				}
				else {
					memcpy(micro_values[i], micro_values[an_individual.prg[i].adr2], micro_params.num_bits_per_dimension * num_dimensions);
					micro_fitness[i] = micro_fitness[an_individual.prg[i].adr2];
				}
				break;

			case MICRO_EA_CROSSOVER:  // Mutation with a fixed mutation probability
				for (int j = 0; j < num_dimensions; j++)
					for (int k = 0; k < micro_params.num_bits_per_dimension; k++) {
						int p = rand() % 2;
						if (p)
							micro_values[i][j * micro_params.num_bits_per_dimension + k] = micro_values[an_individual.prg[i].adr1][j * micro_params.num_bits_per_dimension + k];
						else
							micro_values[i][j * micro_params.num_bits_per_dimension + k] = micro_values[an_individual.prg[i].adr2][j * micro_params.num_bits_per_dimension + k];
					}

				// compute fitness of that micro chromosome
				// transform to base 10
				for (int j = 0; j < num_dimensions; j++)
					x[j] = binary_to_real(micro_values[i], micro_params.num_bits_per_dimension, min_x, max_x);

				micro_fitness[i] = f(x, num_dimensions);// apply f - compute fitness of micro

				break;

			case MICRO_EA_MUTATION:  // Mutation with a fixed mutation probability
				for (int j = 0; j < num_dimensions; j++)
					for (int k = 0; k < micro_params.num_bits_per_dimension; k++) {
						double p = rand() / (double)RAND_MAX;
						if (p < micro_params.mutation_probability)
							micro_values[i][j * micro_params.num_bits_per_dimension + k] = 1 - micro_values[an_individual.prg[i].adr1][j * micro_params.num_bits_per_dimension + k];
						else
							micro_values[i][j * micro_params.num_bits_per_dimension + k] = micro_values[an_individual.prg[i].adr1][j * micro_params.num_bits_per_dimension + k];
					}

				// compute fitness of that micro chromosome
				// transform to base 10
				for (int j = 0; j < num_dimensions; j++)
					x[j] = binary_to_real(micro_values[i], micro_params.num_bits_per_dimension, min_x, max_x);

				micro_fitness[i] = f(x, num_dimensions);// apply f - compute fitness of micro
				break;

			}
		}

		// add to average
		for (int i = 0; i < code_length; i++)
			average_micro_fitness[i] += micro_fitness[i];
	}

	for (int i = 0; i < code_length; i++)
		average_micro_fitness[i] /= (double)micro_params.num_runs;

	an_individual.fitness = average_micro_fitness[0];
	an_individual.best_index = 0;

	for (int i = 1; i < code_length; i++)
	if (an_individual.fitness > average_micro_fitness[i]) {
		an_individual.fitness = average_micro_fitness[i];
		an_individual.best_index = i;
	}


	for (int i = 0; i < code_length; i++) 
		delete[] micro_values[i];

	delete[] micro_values;
	delete[] micro_fitness;
	delete[] average_micro_fitness;
	delete[] x;
}
//---------------------------------------------------------------------------
void start_steady_state_mep(t_meta_gp_parameters &meta_gp_params, t_micro_ea_parameters &micro_params, int num_dimensions, double min_x, double max_x)       // Steady-State 
{
	// a steady state approach:
	// we work with 1 population
	// newly created individuals will replace the worst existing ones (only if they are better).

	// allocate memory
	t_meta_gp_chromosome *population;
	population = new t_meta_gp_chromosome[meta_gp_params.pop_size];
	for (int i = 0; i < meta_gp_params.pop_size; i++)
		allocate_meta_chromosome(population[i], meta_gp_params);

	t_meta_gp_chromosome offspring1, offspring2;
	allocate_meta_chromosome(offspring1, meta_gp_params);
	allocate_meta_chromosome(offspring2, meta_gp_params);


	// initialize
	for (int i = 0; i < meta_gp_params.pop_size; i++) {
		generate_random_meta_chromosome(population[i], meta_gp_params);
		compute_fitness(population[i], meta_gp_params.code_length, micro_params, num_dimensions, min_x, max_x);
	}
	// sort ascendingly by fitness
	qsort((void *)population, meta_gp_params.pop_size, sizeof(population[0]), sort_function);

	printf("generation %d, best fitness = %lf\n", 0, population[0].fitness);

	for (int g = 1; g < meta_gp_params.num_generations; g++) {// for each generation
		for (int k = 0; k < meta_gp_params.pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(population, meta_gp_params.pop_size, 2);
			int r2 = tournament_selection(population, meta_gp_params.pop_size, 2);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < meta_gp_params.crossover_probability)
				one_cut_point_crossover(population[r1], population[r2], meta_gp_params, offspring1, offspring2);
			else {// no crossover so the offspring are a copy of the parents
				copy_individual(offspring1, population[r1], meta_gp_params);
				copy_individual(offspring2, population[r2], meta_gp_params);
			}
			// mutate the result and compute fitness
			mutate_meta_chromosome(offspring1, meta_gp_params);
			compute_fitness(offspring1, meta_gp_params.code_length, micro_params, num_dimensions, min_x, max_x);
			// mutate the other offspring and compute fitness
			mutate_meta_chromosome(offspring2, meta_gp_params);
			compute_fitness(offspring2, meta_gp_params.code_length, micro_params, num_dimensions, min_x, max_x);

			// replace the worst in the population
			if (offspring1.fitness < population[meta_gp_params.pop_size - 1].fitness) {
				copy_individual(population[meta_gp_params.pop_size - 1], offspring1, meta_gp_params);
				qsort((void *)population, meta_gp_params.pop_size, sizeof(population[0]), sort_function);
			}
			if (offspring2.fitness < population[meta_gp_params.pop_size - 1].fitness) {
				copy_individual(population[meta_gp_params.pop_size - 1], offspring2, meta_gp_params);
				qsort((void *)population, meta_gp_params.pop_size, sizeof(population[0]), sort_function);
			}
		}
		printf("generation %d, best fitness = %lf\n", g, population[0].fitness);
	}
	// print best chromosome
	print_meta_chromosome(population[0], meta_gp_params.code_length);

	// free memory
	delete_meta_chromosome(offspring1);
	delete_meta_chromosome(offspring2);

	for (int i = 0; i < meta_gp_params.pop_size; i++)
		delete_meta_chromosome(population[i]);
	delete[] population;
}
//--------------------------------------------------------------------
int main(void)
{
	t_meta_gp_parameters meta_gp_params;

	meta_gp_params.pop_size = 10;						    // the number of individuals in population
	meta_gp_params.code_length = 100;
	meta_gp_params.num_generations = 1000;					// the number of generations
	meta_gp_params.mutation_probability = 0.01;              // mutation probability
	meta_gp_params.crossover_probability = 0.9;             // crossover probability

	t_micro_ea_parameters micro_ea_params;
	micro_ea_params.mutation_probability = 0.01;
	micro_ea_params.num_bits_per_dimension = 30;
	micro_ea_params.num_runs = 30;

	srand(0);
	start_steady_state_mep(meta_gp_params, micro_ea_params, 30, -5, 5);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------