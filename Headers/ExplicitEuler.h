#pragma once

#include <vector>
#include <algorithm>

#include "problem.h"

using namespace std;

/**
 * Runs the Explicit Euler algorithm on a given problem
 *
 * @param problem   The Problem struct to be solved
 * @param N         The number of steps to perform
 * @return          The history of all steps of the method (N x d)
 */
vector<vector<double>> ExplicitEuler(Problem<double>& problem, int N)
{    
    // Step size
    double h = (problem.Tf - problem.t0) / (double)N;
    // Force term
    double* F = new double[problem.d];

    // Prepare the output vector and resize it
    vector<vector<double>> history(N+1);
    for (int i = 0; i < N+1; i++)
        history[i].resize(problem.d);
    
    // Pick the initial data
    for (int i = 0; i < problem.d; i++)
        history[0][i] = problem.u0[i];
    
    // Time stepping
    for (int i = 0; i < N; i++)
    {
        double t = problem.t0 + i * h;
        problem.Force(t, &history[i][0], F, problem.parameters);

        for (int j = 0; j < problem.d; j++)
            history[i+1][j] = history[i][j] + h * F[j];

        // Perform any necessary post processing, clipping or wrapping around of the state variables
        problem.PostProcess(&history[i+1][0]);
    }

    delete[] F;
    return history;
}

