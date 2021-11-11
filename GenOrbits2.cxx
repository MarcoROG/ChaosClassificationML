#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>

#include "Headers/problem.h"
#include "Headers/ProblemLibrary.h"
#include "Headers/Parameters.h"
#include "Headers/MizarHelper.h"
#include "Headers/SplittingMethod.h"
#include "Headers/IO.h"


using namespace std;


/**
 * Computes the tangential flow related to the kinetic part of the Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State of the system
 * @param parameters    Physical parameters for the system
 * @param in            Input tangent vector
 * @param out           Output tangent vector
 */
inline void FluxTangentA(double h, double* state, double* parameters, double* in, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    out[0] = in[0] + w0 * h * in[2];
    out[1] = in[1] + w1 * h * in[3];
    out[2] = in[2];
    out[3] = in[3];    
}

/**
 * Computes the tangential flow related to the potential part of the Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State of the system
 * @param parameters    Physical parameters for the system
 * @param in            Input tangent vector
 * @param out           Output tangent vector
 */
inline void FluxTangentB(double h, double* state, double* parameters, double* in, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    out[0] = in[0];
    out[1] = in[1];

    double F0 = (2.0*q1 + w0) * in[0] + 2.0 * q0 * in[1];
    double F1 = 2.0 * q0 * in[0]  + (- 2.0 * q1 + w1) * in[1];

    out[2] = in[2] - h * F0;
    out[3] = in[3] - h * F1;    
}

/**
 * Computes the tangential flow related to the {{A,B},B} Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State of the system
 * @param parameters    Physical parameters for the system
 * @param in            Input tangent vector
 * @param out           Output tangent vector
 */
inline void FluxTangentCorr(double h, double* state, double* parameters, double* in, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    
    double F0 = (8.0*q0*q0*w1 +8.0*q1*q1*w0 +8.0*q1*w0*w0 +2.0*w0*w0*w0)*in[2] +(8.0*q0*q1*w0 -8.0*q0*q1*w1 -4.0*q0*w0*w0 +4.0*q0*w1*w1)*in[3];
    double F1 = (8.0*q0*q1*w0 -8.0*q0*q1*w1 +4.0*q0*w0*w0 +4.0*q0*w1*w1)*in[2] +(8.0*q0*q0*w0 +8.0*q1*q1*w1 -8.0*q1*w1*w1 +2.0*w1*w1*w1)*in[3];

    out[0] = in[0] + h * F0;
    out[1] = in[1] + h * F1;

    out[2] = in[2];
    out[3] = in[3];    
}


/**
 * Computes the flow related to the kinetic part of the Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State the flow acts upon
 * @param parameters    Physical parameters for the system
 * @param out           Output vector
 */
inline void FluxA(double h, double* state, double* parameters, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    out[0] = state[0] + w0 * h * p0;
    out[1] = state[1] + w1 * h * p1;
    out[2] = state[2];
    out[3] = state[3];    
}

/**
 * Computes the flow related to the potential part of the Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State the flow acts upon
 * @param parameters    Physical parameters for the system
 * @param out           Output vector
 */
inline void FluxB(double h, double* state, double* parameters, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    out[0] = q0;
    out[1] = q1;

    double F0 = - (w0 + 2*q1) * q0;
    double F1 = -w1*q1 -q0*q0 + q1*q1;

    out[2] = state[2] + h * F0;
    out[3] = state[3] + h * F1;
}

/**
 * Computes the flow related to the {{A,B},A} Hamiltonian for HH
 *
 * @param h             Evaluation time
 * @param state         State the flow acts upon
 * @param parameters    Physical parameters for the system
 * @param out           Output vector
 */
inline void FluxCorr(double h, double* state, double* parameters, double* out)
{
    double w0 = parameters[0];
    double w1 = parameters[1];

    double q0 = state[0];
    double q1 = state[1];

    double p0 = state[2];
    double p1 = state[3];

    out[0] = q0;
    out[1] = q1;

    double F0 = -4.0*w1*q0*q0*q0 +(-8.0*w0 +4.0*w1)*q0*q1*q1 +(-8.0*w0*w0 -4.0*w1*w1)*q0*q1 -2.0*w0*w0*w0*q0;
    double F1 = (-8.0*w0+4.0*w1)*q0*q0*q1 +(-4.0*w0*w0-2.0*w1*w1)*q0*q0 -4.0*w1*q1*q1*q1 +6.0*w1*w1*q1*q1 -2.0*w1*w1*w1*q1;

    out[2] = state[2] + h * F0;
    out[3] = state[3] + h * F1;
}

// 1 , 2 , 3 , 4, 5,  6,  7,     8,  9,      10,     11
// w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, id
int main(int argc, char** argv)
{
    // Number of steps to be performed
    double tF = std::stod(argv[8]);

    double Tstep = std::stod(argv[7]);
    unsigned long N = tF / Tstep;

    int solver = std::stoi(argv[10]);
    double Tprint = std::stod(argv[9]);
    unsigned long Nprint = tF / Tprint;

    int SKIP = N / Nprint;

    double t0 = 0;

    double w0 = std::stod(argv[1]);
    double w1 = std::stod(argv[2]);


    // Poincarè plane position
    double q0 = std::stod(argv[3]);

    double* parameters = new double[2];
    parameters[0] = w0;
    parameters[1] = w1;
 
    double* x0 = new double[4];
    double h = tF / (double)N;

    // Determine graphical area
    double Ecrit = std::min(w0*w0*w0/24. + w0*w0*w1/8., w1*w1*w1/6.);
    double E = std::stod(argv[4]);; 
	// Initialize Splitting method
    Corrector* corr = new Corrector();
    vector<double> cs;
    vector<double> ds;
 
    

    SplittingMethod* method = new SplittingMethod();

    if (solver == 0)
    {
        cs = { (5. - sqrt(15.)) / 10. , sqrt(15.) / 10. };
        ds = { 5./18., 4./9. };
        corr->corrector = FluxCorr;
        corr->correctorTang = FluxTangentCorr;
        corr->c = (54. - 13. * sqrt(15.)) / (2. * 648.);
        corr->eps_deg = 2;
        corr->h_deg = 2;
    }
    else
    {
        double sqp = sqrt(490. + 42.*sqrt(105));
        double sqm = sqrt(490. - 42.*sqrt(105));

        cs = { 0.5 - (sqm+sqp)/84., sqm/42., (sqp-sqm)/84.  };     
        ds = { (322.-13.*sqrt(70))/1800., (322.+13*sqrt(70))/1800., 64./225. };
        corr->corrector = FluxCorr;
        corr->correctorTang = FluxTangentCorr;
        corr->c = 0.002270543121419264819434955050039130 / 2.0;
        corr->eps_deg = 2;
        corr->h_deg = 2;   
    }
    method->cs = cs;
    method->ds = ds;
    method->corr = corr;

    Fluxes* fluxes = new Fluxes();
    fluxes->FluxA = FluxA;
    fluxes->FluxB = FluxB;
    fluxes->FluxATang = FluxTangentA;
    fluxes->FluxBTang = FluxTangentB;

    // Setup data    
    x0[0] = q0;
	x0[1] = std::stod(argv[5]);
	x0[3] = std::stod(argv[6]);
    string id = argv[11];

    // Compute p0 such that the energy is fixed
    double q1 = x0[1]; double p1=x0[3];
    double p0square = 2.*E/w0 - (w1/w0)*(p1*p1+q1*q1) - q0*q0 -(2./w0)*q0*q0*q1 + (2./(3.*w0))*q1*q1*q1;
    // If motion is not possible, skip
    if (p0square<0)
        return 0;
    x0[2] = sqrt(p0square);

    Problem<double> prob = ProblemLibrary::QP_HenonHeiles(t0, tF, x0, parameters);

    ///////////////////////////////
    // Actually solve the problem
    ///////////////////////////////
  
    // Initialize result vector
    int ntang = 2;
    double* state = new double[prob.d];
    double* tangFlat = new double[ntang*prob.d];
    double** tangent = new double*[ntang];

    // Note: this makes sense for ntang <= d
    for (int i = 0; i < ntang; i++)
    {
        tangent[i] = new double[prob.d];
        for (int j = 0; j < prob.d; j++)
        {
            tangent[i][j] = (i+j) % prob.d == 0? 1.0 : 0.0;
            tangFlat[i*prob.d + j] = tangent[i][j]; 
        }
    }
      
    for (int i = 0; i < prob.d; i++)
        state[i] = x0[i];

    
    // Open files
    ofstream stateFile("Runs/RUNHIST"+id+".bin", ios::binary | ios::out);
    ofstream tangFile("Runs/RUNTANG"+id+".bin", ios::binary | ios::out);

    // Prepare files
    stateFile << std::fixed << std::setprecision(15); 
    tangFile << std::fixed << std::setprecision(15);  

    int nlines = (N+1) / SKIP;
    int ncolOrb = prob.d;  
    int ncolTang = prob.d * ntang;

    stateFile.write(reinterpret_cast<const char *>(&nlines), sizeof(nlines));
    stateFile.write(reinterpret_cast<const char *>(&ncolOrb), sizeof(ncolOrb));

    tangFile.write(reinterpret_cast<const char *>(&nlines), sizeof(nlines));
    tangFile.write(reinterpret_cast<const char *>(&ncolTang), sizeof(ncolTang));

    
    SaveRowBin(stateFile, state, prob.d);
    SaveRowBin(tangFile, tangFlat, ntang*prob.d);

    for(unsigned long j = 1; j < N+1; j++)
    {
        // Perform a step
        SplittingMethodStepTangent(method, fluxes, state, tangent, ntang, 1.0, h, prob);
        
        // Persist the status
        if (j % SKIP == 0)
        {
            SaveRowBin(stateFile, state, prob.d);

            // Persist the tangent dynamics
            for (int i = 0; i < ntang; i ++)
                std::copy(tangent[i], tangent[i] + prob.d, tangFlat + i*prob.d);

            SaveRowBin(tangFile, tangFlat, prob.d*ntang);

            unsigned int d = prob.d;    
            for (int i = 0; i < ntang; i++)
            {
                double norm = 0.0;
                for (int j = 0; j < d; j++)
                    norm += tangent[i][j] * tangent[i][j];
                norm = sqrt(norm);
                
                for (int j = 0; j < d; j++)
                    tangent[i][j] /= norm;
            }
        }
    }

    delete[] x0;
    return 0;
}
