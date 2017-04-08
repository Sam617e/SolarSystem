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

        cout << "Solar System" << endl;
        Dimension = 3;

        IntegrationPoints = 10000;
        FinalTime = 50.; // units: years

        double TimeStep = FinalTime/((double) IntegrationPoints);
        double x[3],v[3];  // positions and velocities
        //double xJ[3],vJ[3]; // Jupiter positions and velocities (not necessary, but nice if you want to see them)
        // Initial Earth position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        // When adding planets, place them all in a line on the x-axis
        // Always place the Sun last. Idk why, it just works.
        // NOTE (PART C): Try sqrt(GM/r) for speed producing a circular orbit
        // GM = 4*M_PI*M_PI -> v = 2*M_PI/sqrt(r) = 2*M_PI
        // NOTE (PART D): Try sqrt(2GM/r) = sqrt(2)*2*M_PI for escape speed
        // NOTE (PART E): Simply change Jupiter's mass to see how it affects Earth's orbit
        // (mass,x,y,z,vx,vy,vz), Vel ~ 2*M_PI/sqrt(r)
        double r[10]; // contains distance (in AU) from sun for each body, in order of planet distance
        double r[0] = 0; // Sun
        double r[1] = 0.39; // Mercury
        double r[2] = 0.72; // Venus
        double r[3] = 1.0; // Earth
        double r[4] = 1.52; // Mars
        double r[5] = 5.20; // Jupiter
        double r[6] = 9.54; // Saturn
        double r[7] = 19.19; // Uranus
        double r[8] = 30.06; // Neptune
        double r[9] = 39.53 // Pluto

        planet Earth(0.000003,r[3],0,0.0,0,2*M_PI,0);
        planet Jupiter(0.0009546,r[5],0,0,0,2*M_PI/sqrt(5.20),0);
        planet Mercury(0.0000001659,r[1],0,0,0,2*M_PI/sqrt(0.39),0);
        planet Venus(0.00000246,r[2],0,0,0,2*M_PI/sqrt(0.72),0);
        planet Mars(0.0000003318,r[4],0,0,0,2*M_PI/sqrt(1.52),0);
        planet Saturn(0.0002765,r[6],0,0,0,2*M_PI/sqrt(9.54),0);
        planet Uranus(0.00004424,r[7],0,0,0,2*M_PI/sqrt(19.19),0);
        planet Neptune(0.00005178,r[8],0,0,0,2*M_PI/sqrt(30.06),0);
        planet Pluto(0.000000006586,r[9],0,0,0,2*M_PI/sqrt(39.53),0);
        planet Sun(1.,r[]0,0.,0.,0.,0.,0.);

        // NOTE (PART F): Center of Mass calculations for the whole system:
        // int CoM =

        solver binary_vv(5.0);

        // Replicate these lines for each planet you add.
        binary_vv.add(Earth);
        binary_vv.add(Jupiter);
        binary_vv.add(Mercury);
        binary_vv.add(Venus);
        binary_vv.add(Mars);
        binary_vv.add(Saturn);
        binary_vv.add(Uranus);
        binary_vv.add(Neptune);
        binary_vv.add(Pluto);
        binary_vv.add(Sun);

        // Places initial Earth positions and velocities in their arrays
        for(int i = 0; i < Dimension; i++){
            x[i] = Earth.position[i];
            v[i] = Earth.velocity[i];
        }
/*
        for(int i = 0; i < Dimension; i++){
            xJ[i] = Jupiter.position[i];
            vJ[i] = Jupiter.velocity[i];
        }
*/
        cout << "Initial Earth values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);
        cout << endl;
/*
        cout << "Initial Jupiter values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,xJ,vJ,IntegrationPoints);
        cout << endl;
*/
        // Uncomment the second pair of lines to use the Euler algorithm
        cout << "Velocity Verlet results for the system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,1,0.);
        //cout << "Velocity Euler results for the Sun-Earth system:" << endl;
        //binary_vv.VelocityEuler(Dimension,IntegrationPoints,FinalTime,1,0.);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }
/*
        for(int j = 0; j < Dimension;j++){
            xJ[j] = binary_vv.all_planets[1].position[j];
            vJ[j] = binary_vv.all_planets[1].velocity[j];
        }
*/
        PrintFinalValues(Dimension,x,v);
        cout << endl;
        //cout << "Jupiter positions and velocities: " << endl;
        //PrintFinalValues(Dimension,xJ,vJ);

    return 0;
}

/*
// CENTER OF MASS SYSTEM

int main()
{
    int IntegrationPoints;  // No. of integration points
    double FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

        cout << "Earth-Sun binary system" << endl;
        Dimension = 3;

        IntegrationPoints = 10000;
        FinalTime = 50.;
        // time_step = FinalTime/IntegrationPoints

        double TimeStep = FinalTime/((double) IntegrationPoints);
        double x[3],v[3];  // positions and velocities
        double xJ[3],vJ[3]; // Jupiter positions and velocities
        // initial position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        planet planet1(0.000003,0.9950378,0,0.0,0,2*M_PI,0); // Earth: (mass,x,y,z,vx,vy,vz)
        planet planet2(0.0009546,5.195,0,0,0,2*M_PI/sqrt(5.20),0); // Jupiter //vel ~ 2*M_PI/sqrt(5.20) ~ 2.7553590
        planet planet3(1.,-0.0049622,0.,0.,0.,-2*M_PI*(0.000003+0.0009546/sqrt(5.20)),0.); // Sun: (mass,x,y,z,vx,vy,vz)
        //NOTE (PART C): Try sqrt(GM/r) for speed producing a circular orbit
       // GM = 4*M_PI*M_PI -> v = 2*M_PI/sqrt(r) = 2*M_PI
       // NOTE (PART D): Try sqrt(2GM/r) = sqrt(2)*2*M_PI ~ 8.8857659 for escape speed
        solver binary_vv(5.0);
        binary_vv.add(planet1);
        binary_vv.add(planet2);
        binary_vv.add(planet3);

        //cout << "pos: " << x << endl;
        for(int i = 0; i < Dimension; i++){
            x[i] = planet1.position[i];
            v[i] = planet1.velocity[i];
        }

        for(int i = 0; i < Dimension; i++){
            xJ[i] = planet2.position[i];
            vJ[i] = planet2.velocity[i];
        }
        cout << "Initial Earth values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);
        cout << endl;

        cout << "Initial Jupiter values: "<< endl;
        PrintInitialValues(Dimension,TimeStep,FinalTime,xJ,vJ,IntegrationPoints);
        cout << endl;

        cout << "Velocity Verlet results for the Sun-Jupiter system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,1,0.);
        //cout << "Velocity Euler results for the Sun-Earth system:" << endl;
        //binary_vv.VelocityEuler(Dimension,IntegrationPoints,FinalTime,1,0.);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }

        for(int j = 0; j < Dimension;j++){
            xJ[j] = binary_vv.all_planets[1].position[j];
            vJ[j] = binary_vv.all_planets[1].velocity[j];
        }

        PrintFinalValues(Dimension,x,v);
        cout << endl;
        cout << "Jupiter positions and velocities: " << endl;
        PrintFinalValues(Dimension,xJ,vJ);

    return 0;
}
*/
