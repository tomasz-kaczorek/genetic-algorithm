#include "engine.h"

#include <stdlib.h>

#include <QDebug>

Engine::Engine(double (* f)(double const *, size_t), void (* x)(double *, size_t), size_t n, size_t mu, size_t lambda, int maxGenerations, double epsilon, double crossoverProbability, double mutationProbability) :
    m_f(f),
    m_x(x),
    m_n(n),
    m_mu(mu),
    m_lambda(lambda),
    m_generations(0),
    m_maxGenerations(maxGenerations),
    m_epsilon(epsilon),
    m_crossoverProbability(crossoverProbability),
    m_mutationProbability(mutationProbability),
    m_generator(time(nullptr))
{
    std::vector<double> weights(n);
    for(size_t i = 0; i < n; ++i) weights[i] = n - i;
    m_draw = boost::random::discrete_distribution<>(weights.begin(), weights.end());

    m_population = new Chromosome * [mu + lambda];
    for(size_t i = 0; i < mu + lambda; ++i)
    {
        m_population[i] = new Chromosome(n, x);
        rate(m_population[i]);
    }
}

size_t Engine::draw()
{
    return m_draw(m_generator);
}

double Engine::uniform()
{
    return m_uniform(m_generator);
}

double Engine::normal()
{
    return m_normal(m_generator);
}

void Engine::step()
{
    for(size_t i = m_mu; i < m_mu + m_lambda; i += 2)
    {
        duplicate(m_population[i], m_population[draw()]);
        duplicate(m_population[i + 1], m_population[draw()]);
        crossover(m_population[i], m_population[i + 1]);
        mutate(m_population[i]);
        mutate(m_population[i + 1]);
        rate(m_population[i]);
        rate(m_population[i + 1]);
    }
    std::qsort(m_population, m_mu + m_lambda, sizeof(Chromosome *), compare);
}

double Engine::error()
{
    double alleles[m_n];
    double h;
    double t;
    double result;
    result = 0.0;
    memcpy(alleles, m_population[0]->m_alleles, m_n * sizeof(double));
    for(size_t i = 0; i < m_n; ++i)
    {
        t = alleles[i];
        h = sqrt(DBL_EPSILON) * fmax(1.0, fabs(t));
        alleles[i] += h;
        result += pow((m_f(alleles, m_n) - m_f(m_population[0]->m_alleles, m_n)) / h, 2.0);
        alleles[i] = t;
    }
    return sqrt(result);
    return 0.0;
}

void Engine::crossover(Engine::Chromosome * chromosome1, Engine::Chromosome * chromosome2)
{
    if(uniform() < m_crossoverProbability)
    {
        int a = uniform();
        for(size_t i = 0; i < m_n; ++i)
        {
            double allele1 = chromosome1->m_alleles[i];
            double allele2 = chromosome2->m_alleles[i];
            chromosome1->m_alleles[i] = a * allele1 + (1.0 - a) * allele2;
            chromosome2->m_alleles[i] = (1.0 - a) * allele1 + a * allele2;
        }
    }
}

void Engine::mutate(Engine::Chromosome * chromosome)
{
    static double tau1 = 1.0 / sqrt(2.0 * m_n);
    static double tau2 = 1.0 / sqrt(2.0 * sqrt(m_n));

    double tau = tau1 * normal();
    for(size_t i = 0; i < m_n; ++i)
    {
        if(uniform() < m_mutationProbability)
        {
            chromosome->m_alleles[i] = chromosome->m_alleles[i] * chromosome->m_sigma[i] * normal();
            chromosome->m_sigma[i] = chromosome->m_sigma[i] * exp(tau + tau2 * normal());
        }
    }
}

void Engine::rate(Engine::Chromosome * chromosome)
{
    chromosome->m_fitness = m_f(chromosome->m_alleles, m_n);
}

void Engine::duplicate(Engine::Chromosome * destination, Engine::Chromosome * source)
{
    destination->m_fitness = source->m_fitness;
    for(size_t i = 0; i < m_n; ++i)
    {
        destination->m_alleles[i] = source->m_alleles[i];
        destination->m_sigma[i] = source->m_sigma[i];
    }
}

double Engine::nthFitness(size_t n)
{
    return m_population[n]->m_fitness;
}

double * Engine::nthAlleles(size_t n)
{
    return m_population[n]->m_alleles;
}

int Engine::compare(const void * ptr1, const void * ptr2)
{
    Chromosome * chromosome1 = *((Chromosome **) ptr1);
    Chromosome * chromosome2 = *((Chromosome **) ptr2);
    if(chromosome1->m_fitness < chromosome2->m_fitness) return -1;
    if(chromosome1->m_fitness > chromosome2->m_fitness) return 1;
    return 0;
}

Engine::Chromosome::Chromosome(size_t n, void (* alleles)(double *, size_t), void (* sigma)(double *, size_t))
{
    m_alleles = new double[n];
    alleles(m_alleles, n);

    m_sigma = new double[n];
    if(sigma != nullptr) sigma(m_sigma, n);
    else
    {
        for(size_t i = 0; i < n; ++i)
        {
            m_sigma[i] = 1.0;
        }
    }
}

/*

double Engine::Chromosome::error(double (*f)(double * x))
{

}*/
