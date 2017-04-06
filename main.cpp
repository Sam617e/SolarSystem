#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"

using namespace std;

void PrintInitialValues(int, double, double, double *, double *, int);
void PrintFinalValues(int, double *, double *);

void PrintFinalValues(int Dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<Dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<Dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}

void PrintInitialValues(int Dimension,double TimeStep, double FinalTime,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << TimeStep << "; final time = " << FinalTime << "; integration points = " << N << endl;

    cout << "Initial position = ";
    for(int j=0;j<Dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<Dimension;j++) cout << v_initial[j] << " ";
    cout << endl;

}

int main()
{
    int IntegrationPoints;  // No. of integration points
    double FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

        cout << "Earth-Sun binary system" << endl;
        Dimension = 3;

        IntegrationPoints = 10000;
        FinalTime = 500.;
        // time_step = FinalTime/IntegrationPoints

        double TimeStep = FinalTime/((double) IntegrationPoints);
        double x[3],v[3];  // positions and velocities
        // initial position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        planet planet1(0.000003,1.,0.0,0.0,0.0,2*M_PI,0); // Earth: (mass,x,y,z,vx,vy,vz)
        planet planet2(1.,0.,0.,0.,0.,0.,0.);           // Sun: (mass,x,y,z,vx,vy,vz)
        // NOTE (PART C): Try sqrt(GM/r) for speed producing a circular orbit
        // GM = 4*M_PI*M_PI -> v = 2*M_PI/sqrt(r) = 2*M_PI
        solver binary_vv(5.0);
        binary_vv.add(planet1);
        binary_vv.add(planet2);

        //cout << "pos: " << x << endl;
        for(int i = 0; i < Dimension; i++){
            x[i] = planet1.position[i];
            v[i] = planet1.velocity[i];
        }
        //cout << "Initial Earth v: " << planet1.velocity[1] << endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);

        cout << "Velocity Verlet results for the Sun-Earth system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,1,0.);
        //cout << "Velocity Euler results for the Sun-Earth system:" << endl;
        //binary_vv.VelocityEuler(Dimension,IntegrationPoints,FinalTime,1,0.);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }
        PrintFinalValues(Dimension,x,v);
    return 0;
}




