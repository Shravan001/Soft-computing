#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <limits>

// Objective functions
double f1(const std::vector<double>& x) {
    double sum = 0.0;
    for (double value : x) {
        sum += std::pow(value - 1.0 / std::sqrt(x.size()), 2);
    }
    return 1.0 - std::exp(-sum);
}

double f2(const std::vector<double>& x) {
    double sum = 0.0;
    for (double value : x) {
        sum += std::pow(value + 1.0 / std::sqrt(x.size()), 2);
    }
    return 1.0 - std::exp(-sum);
}

// Individual struct
struct Individual {
    std::vector<double> chromosome;
    double fitness1;
    double fitness2;
    double distance;

};

// Function to evaluate objectives for an individual
void evaluate_objectives(Individual& individual) {
    individual.fitness1 = f1(individual.chromosome);
    individual.fitness2 = f2(individual.chromosome);
}

// Random number generator
double random_double(double min, double max) {
    return min + static_cast<double>(rand()) / RAND_MAX * (max - min);
}

// Create a random individual
Individual create_individual(int chromosome_length) {
    Individual individual;
    individual.chromosome.resize(chromosome_length);
    for (int i = 0; i < chromosome_length; ++i) {
        individual.chromosome[i] = random_double(-4.0, 4.0);
    }
    evaluate_objectives(individual);
    return individual;
}

// Crossover operator
Individual crossover(const Individual& parent1, const Individual& parent2) {
    Individual child;
    int size = parent1.chromosome.size();
    child.chromosome.resize(size);
    for (int i = 0; i < size; ++i) {
        double alpha = random_double(0.0, 1.0);
        child.chromosome[i] = alpha * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
    }
    evaluate_objectives(child);
    return child;
}

// Mutation operator
void mutate(Individual& individual, double mutation_rate) {
    for (double& gene : individual.chromosome) {
        if (random_double(0.0, 1.0) < mutation_rate) {
            gene = random_double(-4.0, 4.0);
        }
    }
    evaluate_objectives(individual);
}

// Non-dominated sorting
std::vector<std::vector<Individual> > non_dominated_sort(const std::vector<Individual>& population) {
    std::vector<std::vector<Individual> > fronts(1);
    for (const auto& p : population) {
        bool is_dominated = false;
        for (const auto& q : population) {
            if (&p == &q) continue;
            if (p.fitness1 < q.fitness1 && p.fitness2 < q.fitness2) {
                is_dominated = true;
                break;
            } else if (q.fitness1 < p.fitness1 && q.fitness2 < p.fitness2) {
                is_dominated = false;
                break;
            }
        }
        if (!is_dominated) {
            fronts[0].push_back(p);
        }
    }
    return fronts;
}

// Crowding distance assignment
void crowding_distance_assignment(std::vector<Individual>& front) {
    int size = front.size();
    for (auto& p : front) p.distance = 0.0;
    for (int m = 0; m < 2; ++m) {
        std::sort(front.begin(), front.end(), [&](const Individual& a, const Individual& b) {
            return a.fitness1 < b.fitness1;
        });
        front.front().distance = front.back().distance = std::numeric_limits<double>::infinity();
        double f_max = front.back().fitness1;
        double f_min = front.front().fitness1;
        if (f_max == f_min) continue;
        for (int i = 1; i < size - 1; ++i) {
            front[i].distance += (front[i + 1].fitness1 - front[i - 1].fitness1) / (f_max - f_min);
        }
    }
}

// NSGA-II algorithm
std::vector<Individual> nsga2(int population_size, int chromosome_length, int generations, double crossover_rate, double mutation_rate) {
    std::vector<Individual> population(population_size);
    for (auto& individual : population) {
        individual = create_individual(chromosome_length);
    }
    for (int gen = 0; gen < generations; ++gen) {
        std::vector<Individual> offspring;
        for (int i = 0; i < population_size; i += 2) {
            Individual parent1 = population[i];
            Individual parent2 = population[i + 1];
            if (random_double(0.0, 1.0) < crossover_rate) {
                Individual child1 = crossover(parent1, parent2);
                Individual child2 = crossover(parent2, parent1);
                mutate(child1, mutation_rate);
                mutate(child2, mutation_rate);
                offspring.push_back(child1);
                offspring.push_back(child2);
            }
        }
        std::vector<Individual> combined_population = population;
        combined_population.insert(combined_population.end(), offspring.begin(), offspring.end());
        std::vector<std::vector<Individual> > fronts = non_dominated_sort(combined_population);
        population.clear();
        for (const auto& front : fronts) {
            if (population.size() + front.size() <= population_size) {
                population.insert(population.end(), front.begin(), front.end());
            } else {
                int remaining_space = population_size - population.size();
                std::vector<Individual> sorted_front = front;
                crowding_distance_assignment(sorted_front);
                std::sort(sorted_front.begin(), sorted_front.end(), [](const Individual& a, const Individual& b) {
                    return a.distance > b.distance;
                });
                population.insert(population.end(), sorted_front.begin(), sorted_front.begin() + remaining_space);
                break;
            }
        }
    }
    return population;
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr))); // Seed the random number generator

    int population_size;
    std::cout << "Enter the Population Size: ";
    std::cin >> population_size;

    int chromosome_length;
    std::cout << "\nEnter the Chromosome Length: ";
    std::cin >> chromosome_length;

    int generations;
    std::cout << "\nEnter the number of Generations: ";
    std::cin >> generations;

    double crossover_rate;
    std::cout << "\nEnter the Crossover Rate: ";
    std::cin >> crossover_rate;

    double mutation_rate;
    std::cout << "\nEnter the Mutation Rate: ";
    std::cin >> mutation_rate;

    std::vector<Individual> best_solutions = nsga2(population_size, chromosome_length, generations, crossover_rate, mutation_rate);

    std::cout << "\nBest Solutions:\n";
    for (const auto& solution : best_solutions) {
        std::cout << "Fitness 1: " << solution.fitness1 << " Fitness 2: " << solution.fitness2 << std::endl;
    }

    return 0;
}
