#pragma once

#include <math.h>
#include "problem.h"

double ModTwoPi(double x)
{
    double TwoPi = 2 * M_PI;
    return x-TwoPi * floor(x/TwoPi);
};

/**
 * This static class is used to build IVPs on the go.
 */
class ProblemLibrary
{
public:
    /**
     * Creates an IVP for Henon-Heiles model in Hamiltonian Coordinates
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (q0, q1, p0, p1)
     * @param parameters    Parameters (w0, w1)
     * @return              Problem data structure
     */

    static Problem<double> QP_HenonHeiles(double ti, double tf, double* x0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double w0 = params[0];
            double w1 = params[1];

            res[0] = w0 * u[2];                      // q0 = w0 * p0
            res[1] = w1 * u[3];                      // q1 = w1 * p1
            res[2] = - (w0 + 2*u[1]) * u[0];                     
            res[3] = -w1*u[1] -u[0]*u[0] + u[1]*u[1];
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            return;
        };

        Problem<double> prob(ti, tf, F, x0, 4, parameters, PostProcess);
        return prob;
    }

    /**
     * Creates an IVP for Nonlinear Coupled Oscillators in Hamiltonian Coordinates
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (q0, q1, p0, p1)
     * @param parameters    Parameters (w0, w1)
     * @return              Problem data structure
     */
    static Problem<double> QP_NonlinearCoupledOscillator(double ti, double tf, double* x0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double w0 = params[0];
            double w1 = params[1];

            res[0] = w0 * u[2];                      // q0 = w0 * p0
            res[1] = w1 * u[3];                      // q1 = w1 * p1
            res[2] = 2 * u[1] - (w0 + 2) * u[0];                     
            res[3] = 2 * u[0] - (w1 + 2) * u[1];
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            return;
        };

        Problem<double> prob(ti, tf, F, x0, 4, parameters, PostProcess);
        return prob;
    }

/**
     * Creates an  harmonic oscillator problem in Hamiltonian Coordinates
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (q, p)
     * @param parameters    Parameters (k,m)
     * @return              Problem data structure
     */
    static Problem<double> QP_TwiceHarmonicOscillator(double ti, double tf, double* x0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double k = params[0];
            double m = params[1];

            double w2 = k/m;

            res[0] = 2.*u[1];
            res[1] = - 2. * w2 * u[0];
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            return;
        };
        
        void (*Solution) (double, double*, double*, double*) = [](double t, double* params, double* x0, double* res)
        {
            double q0 = x0[0];
            double p0 = x0[1];            
            double w2 = params[0] / params[1];
            double w = sqrt(w2);

            res[0] = q0 * cos(2.*w*t) + (p0/(2.*w)) * sin(2.*w*t); //q(t)
            res[1] = p0 * cos(2.*w*t) / 2. - w * q0 * sin(2.*w*t); //p(t)
        };
        

        Problem<double> prob(ti, tf, F, x0, 2, parameters, PostProcess, Solution);
        return prob;
    }


    /**
     * Creates an  harmonic oscillator problem in Hamiltonian Coordinates
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (q, p)
     * @param parameters    Parameters (k,m)
     * @return              Problem data structure
     */
    static Problem<double> QP_HarmonicOscillator(double ti, double tf, double* x0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double k = params[0];
            double m = params[1];

            double w = k/m;

            res[0] = u[1];
            res[1] = - w * u[0];
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            return;
        };
        
        void (*Solution) (double, double*, double*, double*) = [](double t, double* params, double* x0, double* res)
        {
            double w = params[0];
            double q0 = x0[0];
            double p0 = x0[1];

            res[0] = q0 * cos(w*t) + (p0/w) * sin(w*t); //q(t)
            res[1] = p0 * cos(w*t) -w * q0 * sin(w*t); //p(t)
        };
        

        Problem<double> prob(ti, tf, F, x0, 2, parameters, PostProcess, Solution);
        return prob;
    }

    /**
     * Creates a coupled oscillators problem in Hamiltonian Coordinates
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (q1, q2, p1, p2)
     * @param parameters    Parameters (k1, k2, k3, m1, m2)
     * @return              Problem data structure
     */
    static Problem<double> QP_CoupledOscillators(double ti, double tf, double* x0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* params)
        {
            double* param = (double*) params;
            double k1 = param[0];
            double k2 = param[1];
            double k3 = param[2];

            double m1 = param[3];
            double m2 = param[4];


            double q1 = u[0];
            double q2 = u[1];
            double p1 = u[2];
            double p2 = u[3];

            res[0] = p1;
            res[1] = p2;
            res[2] = -(k1/m1) * q1 -(k2/m1) * (q1-q2);
            res[3] = -(k3/m2) * q2 -(k2/m2) * (q2-q1);
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            return;
        };
        
        Problem<double> prob(ti, tf, F, x0, 4, parameters, PostProcess);
        return prob;
    }

    /**
     * Creates a simple pendulum problem
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (pos 0-2pi, vel)
     * @param parameters    Parameters (g,l)
     * @return              Problem data structure
     */
    static Problem<double> SimplePendulum(double ti, double tf, double* theta0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double g = params[0];
            double l = params[1];

            res[0] = u[1];
            res[1] = - (g/l) * sin(u[0]);
        };

        
        void (*PostProcess) (double*) = [](double* state)
        {
            state[0] = ModTwoPi(state[0]);
        };

        Problem<double> prob(ti, tf, F, theta0, 2, parameters, PostProcess);
        return prob;
    }

    /**
     * Creates a double pendulum problem
     * 
     * @param ti            Initial time
     * @param tf            Final time
     * @param theta0        Initial data (theta1 0-2pi, theta2 0-2pi, vel1, vel2)
     * @param parameters    Parameters (g,l1, l2, m1, m2)
     * @return              Problem data structure
     */
    static Problem<double> DoublePendulum(double ti, double tf, double* theta0, double* parameters)
    {
        void (*F) (double, double*, double*, void*) = [](double time, double* u, double* res, void* param)
        {
            double* params = (double*) param;
            double g = params[0];
            double l1 = params[1];
            double l2 = params[2];
            double m1 = params[3];
            double m2 = params[4];

            res[0] = u[2];
            res[1] = u[3];

            double theta1 = u[0];
            double theta2 = u[1];

            double deltaTheta1 = u[2];
            double deltaTheta2 = u[3];

            double sin1 = sin(theta1);
            double sin2 = sin(theta2);
            double sinDiff = sin(theta1-theta2);
            double sinTwoDiff = sin( 2.*(theta1-theta2) );
            double cosDiff = cos(theta1-theta2);

            double det = m1 + m2 * sinDiff * sinDiff; // Determinant of A

            // Acceleration of theta_1
            res[2] =            -(m2*l1/2.) * sinTwoDiff    * deltaTheta1 * deltaTheta1;
            res[2] = res[2]     -m2*l2      * sinDiff       * deltaTheta2 * deltaTheta2;
            res[2] = res[2]     -(m1+m2)*g                  * sin1;
            res[2] = res[2]     +m2*g       * cosDiff       * sin2;
            res[2] = res[2] * (l2 / det);
   
            // Acceleration of theta_2
            res[3] =            +(m1+m2)*l1 * sinDiff       * deltaTheta1 * deltaTheta1;
            res[3] = res[3]     +(m2*l2/2.) * sinTwoDiff    * deltaTheta2 * deltaTheta2;
            res[3] = res[3]     +(m1+m2)*g  * cosDiff       * sin1;
            res[3] = res[3]     -(m1+m2)*g                  * sin2;
            res[3] = res[3] * (l1 / det);
        };

        void (*PostProcess) (double*) = [](double* state)
        {
            state[0] = ModTwoPi(state[0]);
            state[1] = ModTwoPi(state[1]);        
        };

        Problem<double> prob(ti, tf, F, theta0, 4, parameters, PostProcess);
        return prob;
    }

private:
    //Hide the constructor to make the class static
    ProblemLibrary() {}
};
