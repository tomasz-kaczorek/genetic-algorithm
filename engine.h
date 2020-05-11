#ifndef ENGINE_H
#define ENGINE_H

#include <boost/random.hpp>

#include <stddef.h>

class Engine
{
private:
    struct Chromosome;

    double (* m_f)(double const *, size_t);
    void (* m_x)(double *, size_t);
    size_t m_n;
    size_t m_mu;
    size_t m_lambda;
    int m_generations;
    int m_maxGenerations;
    double m_epsilon;
    double m_crossoverProbability;
    double m_mutationProbability;
    boost::mt19937 m_generator;
    boost::random::discrete_distribution<> m_draw;
    boost::random::uniform_01<> m_uniform;
    boost::random::normal_distribution<> m_normal;
    Chromosome ** m_population;

public:
    Engine(double (* f)(double const *, size_t), void (* x)(double *, size_t), size_t n, size_t mu = 10, size_t lambda = 20, int maxGenerations = 10000, double epsilon = 0.01, double crossoverProbability = 0.95, double mutationProbability = 0.05);
    size_t draw();
    double uniform();
    double normal();
    void step();
    double error();
    void crossover(Chromosome *, Chromosome *);
    void mutate(Chromosome *);
    void rate(Chromosome *);
    void duplicate(Chromosome *, Chromosome *);
    double nthFitness(size_t n);
    double * nthAlleles(size_t n);
    static int compare(const void *, const void *);
};

struct Engine::Chromosome
{
    double m_fitness;
    double * m_alleles;
    double * m_sigma;
    Chromosome(size_t n, void (* alleles)(double *, size_t), void (* sigma)(double *, size_t) = nullptr);
};

#endif // ENGINE_H
