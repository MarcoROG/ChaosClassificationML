#pragma once

#include <math.h>
#include <vector>
#include <algorithm>

#include "problem.h"

typedef void (*Flux) (double, double*, double*, double*);
typedef void (*FluxTangent) (double, double*, double*, double*, double*);

struct Corrector 
{
    Flux corrector;             /** Poisson bracket for the corrector */
    FluxTangent correctorTang;  /** Poisson bracket for the tangential corrector */
    int h_deg;                  /** Degree for the h in the step */ 
    int eps_deg;                /** Degree for the eps in the step */
    double c;                   /** Constant for the step */
};

struct SplittingMethod
{
    vector<double> cs;
    vector<double> ds;
    Corrector* corr = nullptr;
};

struct Fluxes
{
    Flux FluxA;
    Flux FluxB;
    
    FluxTangent FluxATang = nullptr;
    FluxTangent FluxBTang = nullptr;
};

/**
 * Performs a step of a generic splitting method on the d-dimensional problem H = A + epsilon*B, computing tangential dynamics
 *
 * @param   method      Splitting method to be used for integration
 * @param   fluxes      Fluxes from the hamiltonian splitting
 * @param   state       Current state of the system, used as output too
 * @param   tangent     Tangent vectors
 * @param   ntang       Number of tangent vectors
 * @param   epsilon     Scaling factor for B
 * @param   h           Step size
 * @param   prob        Physical problem
 */
void SplittingMethodStepTangent(SplittingMethod* method, Fluxes* fluxes, double* state, double** tangent, int ntang, double epsilon, double h,
                        Problem<double>& prob)
{    
    Flux FluxA = fluxes->FluxA;
    Flux FluxB = fluxes->FluxB;

    FluxTangent FluxATang = fluxes->FluxATang;
    FluxTangent FluxBTang = fluxes->FluxBTang;

    int steps = method->ds.size();
    int csteps = method->cs.size();
    double corr_const = 0.0;
    Corrector* corr = method->corr;

    // Eventual corrector step
    if (corr != nullptr)
    {
        corr_const = corr->c * pow(epsilon, corr->eps_deg) * pow(h, corr->h_deg);
        for (int i = 0; i < ntang; i++) corr->correctorTang(-h * corr_const, state, prob.parameters, tangent[i], tangent[i]);
        corr->corrector(-h * corr_const, state, prob.parameters, state);
    }

    // c1 d1 ..... c_{steps} d_{steps}
    for(int k = 0; k < steps; k++)
    {
        // Follow A flux for h*c_k
        for (int i = 0; i < ntang; i++) FluxATang(h * method->cs[k], state, prob.parameters, tangent[i], tangent[i]);
        FluxA(h * method->cs[k], state, prob.parameters, state);
        
        // Follow B flux for h*d_k
        for (int i = 0; i < ntang; i++) FluxBTang(epsilon * h * method->ds[k], state, prob.parameters, tangent[i], tangent[i]);
        FluxB(epsilon * h * method->ds[k], state, prob.parameters, state);        
    }   

    // Follow A flux for c_{step}
    for (int i = 0; i < ntang; i++) FluxATang(h * method->cs[csteps-1], state, prob.parameters, tangent[i], tangent[i]);
    FluxA(h * method->cs[csteps-1], state, prob.parameters, state);

    // d{steps-1} c{steps-1} ..... d1 c1
    for(int k = steps-2; k >= 0; k--)
    {
        // Follow B flux for h*d_k
        for (int i = 0; i < ntang; i++) FluxBTang(epsilon * h * method->ds[k], state, prob.parameters, tangent[i], tangent[i]);
        FluxB(epsilon * h * method->ds[k], state, prob.parameters, state);

        // Follow A flux for h*c_k
        for (int i = 0; i < ntang; i++) FluxATang(h * method->cs[k], state, prob.parameters, tangent[i], tangent[i]);
        FluxA(h * method->cs[k], state, prob.parameters, state);        
    }

    // Eventual corrector step
    if (corr != nullptr)
    {
        for (int i = 0; i < ntang; i++) corr->correctorTang(-h * corr_const, state, prob.parameters, tangent[i], tangent[i]);
        corr->corrector(-h * corr_const, state, prob.parameters, state);        
    }

    if (prob.PostProcess != nullptr)
        prob.PostProcess(state);

    /*unsigned int d = prob.d;    
    for (int i = 0; i < ntang; i++)
    {
        double norm = 0.0;
        for (int j = 0; j < d; j++)
            norm += tangent[i][j] * tangent[i][j];
        norm = sqrt(norm);
        
        for (int j = 0; j < d; j++)
            tangent[i][j] /= norm;
    }*/
}

/**
 * Performs a generic splitting method on the d-dimensional problem H = A + epsilon*B with tangent dynamics
 *
 * @param   method      Splitting method to be used for integration
 * @param   fluxes      Fluxes from the hamiltonian splitting
 * @param   epsilon     Scaling factor for B
 * @param   tangentHist History of tangent vectors
 * @param   hist        Result vector
 * @param   prob        Physical problem
 */
void SplittingMethodSolveTangent(SplittingMethod* method, Fluxes* fluxes, double epsilon, vector<vector<double>>& tangentHist, vector<vector<double>>& hist,
                        Problem<double>& prob)
{
    unsigned int d = prob.d;
    int N = hist.size() - 1;
    double h = prob.Tf / (double)N;

    int ntang = tangentHist[0].size() / d; // Number of tangent vectors
    double** tangent = new double*[ntang];

    // Note: this makes sense for ntang <= d
    for (int i = 0; i < ntang; i++)
    {
        tangent[i] = new double[d];
        for (int j = 0; j < d; j++)
        {
            tangent[i][j] = (i+j) % d == 0? 1.0 : 0.0;
            tangentHist[0][i*d + j] = tangent[i][j]; 
        }
    }
        
    double* state = new double[d];

    // Initial data
    for(int i = 0; i < d; i++)
        hist[0][i] = prob.u0[i];

    // Copy initial value in state
    std::copy(prob.u0, prob.u0 + d, state);

    for(int j = 1; j < N+1; j++)
    {
        // Perform a step
        SplittingMethodStepTangent(method, fluxes, state, tangent, ntang, epsilon, h, prob);
        
        // Persist the status
        std::copy(state, state + d, &hist[j][0]);

        // Persist the tangent dynamics
        for (int i = 0; i < ntang; i ++)
            std::copy(tangent[i], tangent[i] + d, &tangentHist[j][i*d]);
    }

    delete[] state;

    for (int i = 0; i < ntang; i++)
        delete[] tangent[i];
    delete[] tangent;
}



/**
 * Performs a step of a generic splitting method on the d-dimensional problem H = A + epsilon*B
 *
 * @param   method      Splitting method to be used for integration
 * @param   fluxes      Fluxes from the hamiltonian splitting
 * @param   state       Current state of the system, used as output too
 * @param   epsilon     Scaling factor for B
 * @param   h           Step size
 * @param   prob        Physical problem
 */
void SplittingMethodStep(SplittingMethod* method, Fluxes* fluxes, double* state, double epsilon, double h,
                        Problem<double>& prob)
{    
    Flux FluxA = fluxes->FluxA;
    Flux FluxB = fluxes->FluxB;

    int steps = method->cs.size();
    double corr_const = 0.0;

    Corrector* corr = method->corr;
    

    // Eventual corrector step
    if (corr != nullptr)
    {
        corr_const = corr->c * pow(epsilon, corr->eps_deg) * pow(h, corr->h_deg);
        corr->corrector(-h * corr_const, state, prob.parameters, state);
    }

    // c1 d1 ..... c_{steps} d_{steps}
    for(int k = 0; k < steps; k++)
    {
        // Follow A flux for h*c_k
        FluxA(h * method->cs[k], state, prob.parameters, state);
        
        // Follow B flux for h*d_k
        FluxB(epsilon * h * method->ds[k], state, prob.parameters, state);
    }   

    // Follow A flux for c_{step}
    FluxA(h * method->cs[steps-1], state, prob.parameters, state);

    // d{steps-1} c{steps-1} ..... d1 c1
    for(int k = steps-2; k >= 0; k--)
    {
        // Follow B flux for h*d_k
        FluxB(epsilon * h * method->ds[k], state, prob.parameters, state);

        // Follow A flux for h*c_k
        FluxA(h * method->cs[k], state, prob.parameters, state);
    }

    // Eventual corrector step
    if (corr != nullptr)
    {
        corr->corrector(-h * corr_const, state, prob.parameters, state);
    }

    if (prob.PostProcess != nullptr)
        prob.PostProcess(state);
}


/**
 * Performs a generic splitting method on the d-dimensional problem H = A + epsilon*B
 *
 * @param   method      Splitting method to be used for integration
 * @param   fluxes      Fluxes from the hamiltonian splitting
 * @param   epsilon     Scaling factor for B
 * @param   hist        Result vector
 * @param   prob        Physical problem
 */
void SplittingMethodSolve(SplittingMethod* method, Fluxes* fluxes, double epsilon, vector<vector<double>>& hist,
                        Problem<double>& prob)
{
    unsigned int d = prob.d;
    int N = hist.size() - 1;
    double h = prob.Tf / (double)N;
        
    double* state = new double[d];

    // Initial data
    for(int i = 0; i < d; i++)
        hist[0][i] = prob.u0[i];

    std::copy(prob.u0, prob.u0 + d, state);

    for(int j = 1; j < N+1; j++)
    {
        SplittingMethodStep(method, fluxes, state, epsilon, h, prob);        
        // Persist the status
        std::copy(state, state + d, &hist[j][0]);
    }

    delete[] state;
}

