#ifndef PROBLEM_H_INCLUDED
#define PROBLEM_H_INCLUDED

template <typename T>
/**
 * A struct representing an IVP
 */
struct Problem
{
    void (*Force) (T, T*, T*, void*); /*!< The function calculating the force on members. Accepts: the time, a status vector, a result vector, additional parameters */
    void (*PostProcess) (double*); /*!< Post process and validate the result at each step */
    void (*ExactSolution) (double, double*, double*, double*); /*!< Exact solution if available: Arguments are (time, parameters, initial data, result) */
    T t0; /*!< Initial time of the IVP */
    T Tf; /*!< Final time of the IVP */
    T* u0; /*!< Initial state of the IVP */
    double* parameters; /*!< Parameters and constants for the problem */
    unsigned int d; /*!< Dimensionality of the problem */


    /**
     * Build an IVP with the required input data
     * @param t0_p      The starting time
     * @param tf_p      The ending time
     * @param Force_p   The force acting on the system
     * @param u0_p      The starting condition of the system
     * @param d_p       The problem dimension after order reduction
     * @param p         Problem parameters and costants
     * @param post_proc Function to validate and post process each step
     */
    Problem (T t0_p, T tf_p, void (*Force_p) (T, T*, T*, void*), T* u0_p, unsigned int d_p, double* p, void (*post_proc) (double*)) 
        : Force(Force_p), t0(t0_p), Tf(tf_p), u0(u0_p), d(d_p), parameters(p), PostProcess(post_proc), ExactSolution(nullptr){};

    /**
     * Build an IVP with the required input data and a given solution
     * @param t0_p      The starting time
     * @param tf_p      The ending time
     * @param Force_p   The force acting on the system
     * @param u0_p      The starting condition of the system
     * @param d_p       The problem dimension after order reduction
     * @param p         Problem parameters and costants
     * @param post_proc Function to validate and post process each step
     * @param sol       Exact solution of the problem
     */
    Problem (T t0_p, T tf_p, void (*Force_p) (T, T*, T*, void*), T* u0_p, unsigned int d_p, double* p, void (*post_proc) (double*), void (*sol) (double, double*, double*, double*)) 
        : Force(Force_p), t0(t0_p), Tf(tf_p), u0(u0_p), d(d_p), parameters(p), PostProcess(post_proc), ExactSolution(sol){};
};

#endif // PROBLEM_H_INCLUDED
