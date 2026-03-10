#include "../include/funciones_1d.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <omp.h>
#include <tuple>

void hello_world_funciones(){
    cout<<"Hello, world! (funciones)"<<endl;
}

vector<double> elem_1d_2n::Ns(const vector<double> rho) const{
    return {{(1.0-rho[0])/2.0, (1.0+rho[0])/2.0}};
}

vector<vector<double> > elem_1d_2n::DNs(const vector<double> rho, bool transposicion) const{
    if(transposicion) return vector<vector<double> > {{-0.5, 0.5}};
    else return vector<vector<double> > {{-0.5}, {0.5}};
}

vector<double> elem_1d_3n::Ns(const vector<double> rho) const{
    return {{(rho[0]-pow(rho[0], 2.0))/2.0, 1-pow(rho[0], 2.0), (rho[0]+pow(rho[0], 2.0))/2.0}};
}

vector<vector<double> > elem_1d_3n::DNs(const vector<double> rho, bool transposicion) const{
    if(transposicion) return vector<vector<double> > {{(2*rho[0]-1.0)/2.0, -2.0*rho[0], (2*rho[0]+1.0)/2.0}};
    else return vector<vector<double> > {{(2*rho[0]-1.0)/2.0}, {-2.0*rho[0]}, {(2*rho[0]+1.0)/2.0}};
}

/*
LEGACY
vector<vector<double> > mat_calor(int n, double l, double k){
    double deltax=l/n;
    vector<vector<double> > vect(n, vector<double>(n, 0.0));

    for(int i=0;i<n;i++){
        if(i>0) vect[i-1][i]=vect[i][i-1]=-1.0*(k/deltax);
        vect[i][i]=2.0*(k/deltax);
    }
    vect[n-1][n-1]=1*(k/deltax);

    return vect;
}

vector<double> generar_vector_eccalor(int n, double l, double k, double q, double phi0){
    vector<double> vect(n);
    double deltax=l/n, element=q*deltax;

    for(int i=1;i<vect.size();i++) vect[i]=element;

    vect[0]=(element)+(phi0*k)/deltax;
    vect[vect.size()-1]=element/2.0;
    
    return vect;
}
*/