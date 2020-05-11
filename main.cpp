#include "engine.h"

#include <QDebug>

double f(double const * X, size_t n)
{
    double result = 0.0;
    for(size_t i = 2; i < n; ++i)
    {
        result += 100 * (pow(X[i], 2) + pow(X[i-1], 2)) + pow(X[i-2], 2);
    }
    return result;
}

void x(double * X, size_t n)
{
    for(size_t i = 0; i < n; ++i)
    {
        X[i] = 3.0;
    }
}


int main(int argc, char *argv[])
{
    Engine engine(f, x, 10);
    for(int i = 0; i < 1000; ++i)
    {
        engine.step();
        qDebug("Function evaluatiob: %f\nError: %f\nAlleles: ", engine.nthFitness(0), engine.error());
        for(size_t j = 0; j < 10; ++j)
        {
            qDebug("%f ", engine.nthAlleles(0)[j]);
        }
        qDebug("\n");
    }
    return 0;
}
