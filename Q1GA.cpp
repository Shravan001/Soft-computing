#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include<limits>
#include <algorithm>

using namespace std;

// Define the objective functions f1 and f2
double f1(double x, double y) {
    return 2 + pow(x - 2, 2) + pow(y - 1, 2);
}

double f2(double x, double y) {
    return 9 * x - pow(y - 1, 2);
}

// Combined objective function
double combined_objective(double x, double y) {
    return 0.5 * f1(x, y) + 0.5 * f2(x, y);
}

// Define chromosome representation
vector<double> create_chromosome(int length) {
    vector<double> chromosome(length);
    for (int i = 0; i < length; ++i) {
        chromosome[i] = -20 + static_cast<double>(rand()) / RAND_MAX * 40; // Random value between -20 and 20
    }
    return chromosome;
}

// Fitness function
double fitness(const vector<double>& chromosome) {
    double x = chromosome[0];
    double y = chromosome[1];
    return -combined_objective(x, y);
}

// Tournament selection
vector<double> tournament_selection(const vector<vector<double> >& population, const vector<double>& fitness_values, int tournament_size) {
    vector<double> best_chromosome;
    double best_fitness = -numeric_limits<double>::infinity();
    for (int i = 0; i < tournament_size; ++i) {
        int index = rand() % population.size();
        if (fitness_values[index] > best_fitness) {
            best_fitness = fitness_values[index];
            best_chromosome = population[index];
        }
    }
    return best_chromosome;
}

// Arithmetic crossover
vector<double> crossover(const vector<double>& parent1, const vector<double>& parent2) {
    vector<double> child(parent1.size());
    for (int i = 0; i < parent1.size(); ++i) {
        double alpha = static_cast<double>(rand()) / RAND_MAX;
        child[i] = alpha * parent1[i] + (1.0 - alpha) * parent2[i];
    }
    return child;
}

// Uniform mutation
void mutate(vector<double>& chromosome, double mutation_rate) {
    for (int i = 0; i < chromosome.size(); ++i) {
        if (static_cast<double>(rand()) / RAND_MAX < mutation_rate) {
            chromosome[i] = -20 + static_cast<double>(rand()) / RAND_MAX * 40; // Random value between -20 and 20
        }
    }
}

// Genetic Algorithm
vector<double> genetic_algorithm(int population_size, int chromosome_length, int tournament_size, double mutation_rate, int generations) {
    vector<vector<double> > population(population_size);
    for (int i = 0; i < population_size; ++i) {
        population[i] = create_chromosome(chromosome_length);
    }

    for (int gen = 0; gen < generations; ++gen) {
        // Evaluate fitness
        vector<double> fitness_values(population_size);
        for (int i = 0; i < population_size; ++i) {
            fitness_values[i] = fitness(population[i]);
        }

        // Select parents
        vector<vector<double> > parents(population_size);
        for (int i = 0; i < population_size; ++i) {
            parents[i] = tournament_selection(population, fitness_values, tournament_size);
        }

        // Crossover
        vector<vector<double> > offspring;
        for (int i = 0; i < population_size; i += 2) {
            offspring.push_back(crossover(parents[i], parents[i + 1]));
        }

        // Mutation
        for (auto& child : offspring) {
            mutate(child, mutation_rate);
        }

        // Replace population
        population = offspring;
    }

    // Find the best solution
    double best_fitness = -numeric_limits<double>::infinity();
    vector<double> best_chromosome;
    for (const auto& chromosome : population) {
        double current_fitness = fitness(chromosome);
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            best_chromosome = chromosome;
        }
    }
    return best_chromosome;
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr))); // Seed the random number generator

    int population_size;
    cout << "Enter the Population Size: ";
    cin >> population_size;

    int chromosome_length = 2; // Two parameters x and y
    int tournament_size;
    cout << "\nEnter the Tournament Size: ";
    cin >> tournament_size;

    double mutation_rate;
    cout << "\nEnter the Mutation Rate: ";
    cin >> mutation_rate;

    int generations;
    cout << "\nEnter the number of Generations: ";
    cin >> generations;

    vector<double> best_solution = genetic_algorithm(population_size, chromosome_length, tournament_size, mutation_rate, generations);
    double best_fitness = fitness(best_solution);

    cout << "\nBest Solution: ";
    for (double value : best_solution) {
        cout << value << " ";
    }
    cout << endl;
    cout << "\nBest Fitness: " << best_fitness << endl;

    return 0;
}
