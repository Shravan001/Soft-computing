#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

using namespace std;

// Structure to represent an individual in the population
struct Individual {
    double x; // x-coordinate
    double y; // y-coordinate
    double fitness1; // Fitness value for the first objective function
    double fitness2; // Fitness value for the second objective function
    double crowding_distance; // Crowding distance
    int rank; // Non-dominance rank
};

// Objective functions
double f1(double x, double y) {
    return 2 + pow(x - 2, 2) + pow(y - 1, 2);
}

double f2(double x, double y) {
    return 9 * x - pow(y - 1, 2);
}

// Constraint functions
bool g1(double x, double y) {
    return (pow(x, 2) + pow(y, 2) <= 225);
}

bool g2(double x, double y) {
    return (x - 3 * y + 10 <= 0);
}

// Generate a random number between min and max
double random_double(double min, double max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(min, max);
    return dis(gen);
}

// Initialize the population with random individuals within the search domain
vector<Individual> initialize_population(int pop_size) {
    vector<Individual> population(pop_size);
    for (int i = 0; i < pop_size; ++i) {
        population[i].x = random_double(-20, 20);
        population[i].y = random_double(-20, 20);
        population[i].fitness1 = f1(population[i].x, population[i].y);
        population[i].fitness2 = f2(population[i].x, population[i].y);
    }
    return population;
}

// Perform tournament selection to select individuals for mating
Individual tournament_selection(const vector<Individual>& population) {
    int idx1 = rand() % population.size();
    int idx2 = rand() % population.size();
    return (population[idx1].fitness1 < population[idx2].fitness1 &&
            population[idx1].fitness2 < population[idx2].fitness2) ? population[idx1] : population[idx2];
}

// Perform crossover between two parents to generate a child
Individual crossover(const Individual& parent1, const Individual& parent2, double crossover_prob) {
    Individual child;
    if (random_double(0, 1) < crossover_prob) {
        child.x = (parent1.x + parent2.x) / 2.0;
        child.y = (parent1.y + parent2.y) / 2.0;
    } else {
        child = (random_double(0, 1) < 0.5) ? parent1 : parent2;
    }
    return child;
}

// Perform mutation on an individual
void mutate(Individual& individual, double mutation_prob) {
    if (random_double(0, 1) < mutation_prob) {
        individual.x += random_double(-1, 1);
        individual.y += random_double(-1, 1);
    }
}

// Assign crowding distance to individuals in a front
void assign_crowding_distance(vector<Individual>& front) {
    int n = front.size();
    for (int i = 0; i < n; ++i) {
        front[i].crowding_distance = 0.0;
    }
    for (int m = 0; m < 2; ++m) {
        sort(front.begin(), front.end(), [m](const Individual& a, const Individual& b) {
            return (m == 0) ? (a.fitness1 < b.fitness1) : (a.fitness2 < b.fitness2);
        });
        double f_max = (m == 0) ? front.back().fitness1 : front.back().fitness2;
        double f_min = (m == 0) ? front.front().fitness1 : front.front().fitness2;
        if (f_max == f_min) continue;
        front.front().crowding_distance = front.back().crowding_distance = numeric_limits<double>::infinity();
        for (int i = 1; i < n - 1; ++i) {
            front[i].crowding_distance += (m == 0) ? (front[i + 1].fitness1 - front[i - 1].fitness1) / (f_max - f_min)
                                                    : (front[i + 1].fitness2 - front[i - 1].fitness2) / (f_max - f_min);
        }
    }
}

// Main NSGA-II algorithm
vector<Individual> nsga2(int pop_size, int num_generations, double crossover_prob, double mutation_prob) {
    vector<Individual> population = initialize_population(pop_size);
    for (int generation = 0; generation < num_generations; ++generation) {
        vector<Individual> offspring;
        for (int i = 0; i < pop_size; ++i) {
            Individual parent1 = tournament_selection(population);
            Individual parent2 = tournament_selection(population);
            Individual child = crossover(parent1, parent2, crossover_prob);
            mutate(child, mutation_prob);
            child.fitness1 = f1(child.x, child.y);
            child.fitness2 = f2(child.x, child.y);

            // Ensure the child satisfies the constraints
            if (!g1(child.x, child.y) || !g2(child.x, child.y)) {
                // If the child violates constraints, assign a high fitness to discourage its selection
                child.fitness1 = child.fitness2 = numeric_limits<double>::max();
            }
            offspring.push_back(child);
        }
        population.insert(population.end(), offspring.begin(), offspring.end());
        // Non-dominated sorting
        vector<vector<Individual> > fronts;
        vector<Individual> front;
        for (auto& individual : population) {
            individual.rank = 0;
            for (auto& other : population) {
                if (&individual == &other) continue;
                if (individual.fitness1 > other.fitness1 && individual.fitness2 > other.fitness2) {
                    individual.rank++;
                }
            }
            if (individual.rank == 0) {
                front.push_back(individual);
            }
        }
        assign_crowding_distance(front);
        fronts.push_back(front);
        population.clear();
        for (auto& f : fronts) {
            population.insert(population.end(), f.begin(), f.end());
            if (population.size() >= pop_size) {
                break;
            }
        }
        if (population.size() > pop_size) {
            population.resize(pop_size);
        }
    }
    return population;
}

// Main function
int main() {
    // Input parameters from the user
    int pop_size, num_generations;
    double crossover_prob, mutation_prob;

    cout << "Enter population size: ";
    cin >> pop_size;

    cout << "Enter number of generations: ";
    cin >> num_generations;

    cout << "Enter crossover probability (0-1): ";
    cin >> crossover_prob;

    cout << "Enter mutation probability (0-1): ";
    cin >> mutation_prob;

    // Run NSGA-II algorithm
    vector<Individual> solution_set = nsga2(pop_size, num_generations, crossover_prob, mutation_prob);

    // Print the final solution set
    cout << "Final Solution Set:\n";
    for (const auto& individual : solution_set) {
        cout << "x: " << individual.x << " y: " << individual.y << " f1: " << individual.fitness1 << " f2: " << individual.fitness2 << endl;
    }
    return 0;
}
