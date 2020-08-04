/*
Euler forward method vs Euler backward method vs Crank-Nichlson method vs RK2 method.
Author: Tomas Coronado Gonzalez
Date: July 2020

IVP:   {ODE: y'(t) = y(t)*cos(t)
       {     y(0) = 1

The explicit methods have a numerical inestability (solution -> infinity, 
with t -> infinity) if the step size is big. However, compared with 
implicit methods, they are much easier to implement.

In order to obtain good solutions first-order methods require finer step 
sizes than second-order methods.
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double f (double w_fun, double t_fun_w);         // dy/dt=f(y,t) ODE expressed in its normal form
double f_sol (double t_fun_sol);                 // Exact solution
double Newton_Raphson_w_impl (double w0,         // Newton-Raphson algorithm for solving the implicit ecuation in Euler forward method
                              double w_previous, 
                              double t_now, 
                              double h_interval); 
double Newton_Raphson_w_cn (double w0,           // Newton-Raphson algorithm for solving the implicit ecuation in C-N method
                            double w_previous,
                            double t_previous, 
                            double t_now, 
                            double h_interval);                                                  

int main()
{
    FILE *pf;
    
    int n=1000;                                  // Number of steps of the interval
    int j,k;
    
    double a = 0;                                // Begining of the interval
    double b = 50;                               // End of the interval             
    double h = (b-a)/n;                          // Lengh of a subinterval
    
    double t[n+1];                               // Time array
    t[0] = a;
    double y[n+1];                               // Exact solution array
    y[0] = 1;

    // The next variables are the aproximate solution for each method:
    double w_expl[n+1];
    w_expl[0] = y[0];
    double w_impl[n+1];
    w_impl[0] = y[0];
    double w_cn[n+1];
    w_cn[0] = y[0];
    double w_RK2[n+1], k1, k2;
    w_RK2[0] = y[0];

    for (j=1;j<n+1;j++)
    {
        t[j] = 0;
        y[j] = 0;
        w_expl[j] = 0;
        w_impl[j] = 0;
        w_cn[j]   = 0;
        w_RK2[j]  = 0;
    }
    
    for (k=1;k<n+1;k++)              
    {
        t[k] = t[k-1] + h; 

        // Euler forward method:                        
        w_expl[k] = w_expl[k-1] + h*f(w_expl[k-1],t[k-1]);
        
        // Euler backward method:                // We need to solve w_impl[k]=w_impl[k-1]+h*f(w_impl[k],t[k]), for w_impl[k]
        w_impl[k] = Newton_Raphson_w_impl(w_impl[k-1], w_impl[k-1], t[k], h);
        
        // Runge-Kutta with two stages (RK2):
        k1 = f(w_RK2[k-1],t[k-1]);
        k2 = f(w_RK2[k-1]+h*k1,t[k]);
        w_RK2[k]=w_RK2[k-1] + h*(0.5*k1+0.5*k2);

        // Crank Nicolson method:                // We need to solve w_cn[k]=w_cn[k-1]+0.5*h*f(w_cn[k-1],t[k-1])+0.5*h*(w_cn[k],t[k]), for w_cn[k]
        w_cn[k] = Newton_Raphson_w_cn(w_cn[k-1], w_cn[k-1], t[k-1], t[k], h);                                                                
    
        // Exact solution:
        y[k] = f_sol(t[k]);
    }
    
    printf("Creating file...\n"); 

    // Open file:
    //pf = fopen("output_P1.txt", "w");
    pf = fopen("output_P1.csv", "w");
    if (pf == NULL)
    {
        printf("Error: the file cannot be opened\n");
        exit(1);
    }
    else
    { 
        // Writting on the file:
        fprintf(pf, "time,w_expl,w_impl,w_RK2,w_cn,exact_solution\n");
        for (j=0;j<n+1;j++)
        {
            // fprintf(pf, "%f %f %f %f %f\n", t[j], w_expl[j], w_impl[j], w_RK2[j], y[j]);
            fprintf(pf, "%f,%f,%f,%f,%f,%f\n", t[j], w_expl[j], w_impl[j], w_RK2[j], w_cn[j], y[j]);
        }
        // Close file:
        fclose(pf);
        printf("File with the solution (created and closed)\n");
    }
    if (ferror(pf))
    {
        printf("Error when creating the file\n");
        clearerr(pf);
    }
    
    return 0;
}

double f (double w_fun, double t_fun_w) {
    double funcion_w;
    funcion_w = w_fun*cos(t_fun_w);

    return(funcion_w);
}

double f_sol (double t_fun_sol) {
    double funcion_sol;
    funcion_sol = exp(sin(t_fun_sol));

    return(funcion_sol);
}

double Newton_Raphson_w_impl (double w0, double w_previous, double t_now, double h_interval) {
    // w0 is the first iteraction. And in this code we are using w_impl[k-1] as the first iteraction to calculate w_impl[k]
    int itr = 1, q = 0, maxitr = 100;
    double h_nr, w1, sol = 0.0, allerr = 1e-3; // allowed error
    do
    {
        // h_nr = f(w0)/df(w0);    // W = w_previous + h_interval * W * cos(t_now);  // W = w_previous + h_interval * f(W,t_now);
                                   // function -> f : W - (w_previous + h_interval * W * cos(t_now))
                                   // derivate -> df : 1 - (h_interval * cos(t_now))
        h_nr = (w0 - (w_previous + h_interval * w0 * cos(t_now))) / (1 - (h_interval * cos(t_now)));
        w1 = w0 - h_nr;
        if (fabs(h_nr) < allerr && q==0)
        {
            sol = w1;
            q++;
        }
        w0=w1;

        itr++;
    } while (itr<=maxitr && q==0);
    if (q==0) printf(" The required solution for EB does not converge or iterations are insufficient\n");

    return sol;
}

double Newton_Raphson_w_cn (double w0, double w_previous, double t_previous, double t_now, double h_interval) {
    // w0 is the first iteraction. And in this code we are using w_impl[k-1] as the first iteraction to calculate w_impl[k]
    int itr = 1, q = 0, maxitr = 100;
    double h_nr, w1, sol = 0.0, allerr = 1e-3; // allowed error
    do
    {
        // h_nr = f(w0)/df(w0);    // W = w_previous + 0.5 * h_interval * w_previous * cos(t_previous) + 0.5 * h_interval * W * cos(t_now);  // W - (w_previous + 0.5 * h * f(w_previous,t_previous) + 0.5 * h * f(W,t_now));  
                                   // function -> f : W - (w_previous + 0.5 * h_interval * w_previous * cos(t_previous) + 0.5 * h_interval * W * cos(t_now));
                                   // derivate -> df : 1 - (0.5 * h_interval * cos(t_now))
        h_nr = (w0 - (w_previous + 0.5 * h_interval * w_previous * cos(t_previous) + 0.5 * h_interval * w0 * cos(t_now))) / (1 - (0.5 * h_interval * cos(t_now)));
        w1 = w0 - h_nr;
        if (fabs(h_nr) < allerr && q==0)
        {
            sol = w1;
            q++;
        }
        w0=w1;

        itr++;
    } while (itr<=maxitr && q==0);
    if (q==0) printf(" The required solution for CN does not converge or iterations are insufficient\n");

    return sol;
}





