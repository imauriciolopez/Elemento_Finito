#include "funciones_fem.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>

vector<double> factorizacion_choleski_llt(vector<vector<double> > arri, vector<double> vect) {
    int i, j, k, n, n_elements=arri.size();
    double suma;

    vector<vector<double> > mat(n_elements, vector<double>(n_elements, 0.0));
    vector<double> be(n_elements, 0.0);
    vector<double> salida(n_elements, 0.0);
    
    for(i=0;i<n_elements;i++){
        for(j=0;j<=i;j++){
            suma=0.0;
            for(k=0;k<j;k++){
                suma+=mat[i][k]*mat[j][k];
            }
            if(i==j){
                if(arri[i][i]-suma<0){
                    cout<<"EEEERRRRROR"<<endl;
                }
                mat[i][j]=sqrt(arri[i][i]-suma);
            }
            else{
                mat[i][j]=(arri[i][j]-suma)/mat[j][j];
            }
        }
    }

    for(i=0;i<n_elements;i++){
        suma=0.0;
        for(j=0;j<i;j++){
            suma+=mat[i][j]*be[j];
        }
        be[i]=(vect[i]-suma)/mat[i][i];
    }

    for(i=n_elements-1;i>=0;i--){
        suma=0.0;
        for(j=i+1;j<n_elements;j++){
            suma+=mat[j][i]*salida[j];
        }
        salida[i]=(be[i]-suma)/mat[i][i];
    }

    return salida;
}

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

vector<vector<double> > matriz_rigidez_elemental(){
    return {};
}