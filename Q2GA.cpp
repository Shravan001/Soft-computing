#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>

// Define the objective functions f1 and f2
double f1(const std::vector<double>& x) {
    double sum = 0.0;
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        sum += std::pow(x[i] - 1.0 / std::sqrt(n), 2);
    }
    return 1.0 - std::exp(-sum);
}

double f2(const std::vector<double>& x) {
    double sum = 0.0;
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        sum += std::pow(x[i] + 1.0 / std::sqrt(n), 2);
    }
    return 1.0 - std::exp(-sum);
}

// Combined objective function
double combined_objective(const std::vector<double>& x) {
    return 0.5 * f1(x) + 0.5 * f2(x);
}

// Define chromosome representation
std::vector<double> create_chromosome(int length) {
    std::vector<double> chromosome(length);
    for (int i = 0; i < length; ++i) {
        chromosome[i] = -4.0 + static_cast<double>(rand()) / RAND_MAX * 8.0; // Random value between -4 and 4
    }
    return chromosome;
}

// Fitness function
double fitness(const std::vector<double>& chromosome) {
    return -combined_objective(chromosome);
}

// Tournament selection
std::vector<double> tournament_selection(const std::vector<std::vector<double> >& population, const std::vector<double>& fitness_values, int tournament_size) {
    std::vector<double> best_chromosome;
    double best_fitness = -std::numeric_limits<double>::infinity();
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
std::vector<double> crossover(const std::vector<double>& parent1, const std::vector<double>& parent2) {
    std::vector<double> child(parent1.size());
    for (int i = 0; i < parent1.size(); ++i) {
        double alpha = static_cast<double>(rand()) / RAND_MAX;
        child[i] = alpha * parent1[i] + (1.0 - alpha) * parent2[i];
    }
    return child;
}

// Uniform mutation
void mutate(std::vector<double>& chromosome, double mutation_rate) {
    for (int i = 0; i < chromosome.size(); ++i) {
        if (static_cast<double>(rand()) / RAND_MAX < mutation_rate) {
            chromosome[i] = -4.0 + static_cast<double>(rand()) / RAND_MAX * 8.0; // Random value between -4 and 4
        }
    }
}

// Genetic Algorithm
std::vector<double> genetic_algorithm(int population_size, int chromosome_length, int tournament_size, double mutation_rate, int generations) {
    std::vector<std::vector<double> > population(population_size);
    for (int i = 0; i < population_size; ++i) {
        population[i] = create_chromosome(chromosome_length);
    }

    for (int gen = 0; gen < generations; ++gen) {
        // Evaluate fitness
        std::vector<double> fitness_values(population_size);
        for (int i = 0; i < population_size; ++i) {
            fitness_values[i] = fitness(population[i]);
        }

        // Select parents
        std::vector<std::vector<double> > parents(population_size);
        for (int i = 0; i < population_size; ++i) {
            parents[i] = tournament_selection(population, fitness_values, tournament_size);
        }

        // Crossover
        std::vector<std::vector<double> > offspring;
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
    double best_fitness = -std::numeric_limits<double>::infinity();
    std::vector<double> best_chromosome;
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
    std::cout<<"Enter the Population Size: ";
    std::cin>>population_size;

    int chromosome_length;
    std::cout<<"\nEnter the Chrmosome Length: ";
    std::cin>>chromosome_length;

    int tournament_size;
    std::cout<<"\nEnter the Tournament Size: ";
    std::cin>>tournament_size;

    double mutation_rate;
    std::cout<<"\nEnter the Mutation Rate: ";
    std::cin>>mutation_rate;

    int generations;
    std::cout<<"\nEnter the number of Generations: ";
    std::cin>>generations;

    std::vector<double> best_solution = genetic_algorithm(population_size, chromosome_length, tournament_size, mutation_rate, generations);
    double best_fitness = fitness(best_solution);

    std::cout << "\nBest Solution: ";
    for (double value : best_solution) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    std::cout << "\nBest Fitness: " << best_fitness << std::endl;

    return 0;
}