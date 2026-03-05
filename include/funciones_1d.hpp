#ifndef FUNCIONES_HPP
#define FUNCIONES_HPP

#include "../include/utils.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <tuple>

using namespace std;

void hello_world_funciones();

class elem_1d_2n:public elemento{
public:
    //                    ┌1 dim
    //                    |  ┌2 nodos
    elem_1d_2n(vector<int> indices, vector<struct nodo*> nodos_, vector<vector<double> > D, double cond_neumann_=0.0):elemento(1, 2, "lineal"){
        nodos=nodos_;
        nodos_posiciones=vector<vector<double> >(2, vector<double>(1, 0.0));
        for(int i=0;i<2;i++){
            nodos_posiciones[i]=nodos[i]->posicion;
        }
        material.D=D;
        material.cond_neuman=cond_neumann_;
    }

    vector<double> Ns(const vector<double> rho) const override;

    vector<vector<double> > DNs(const vector<double> rho, bool transposicion=true) const override;

    vector<vector<double> > matriz_rigidez_elemental() const override;

    vector<double> vector_fuerza_elemental() const override;
};



/*
LEGACY
vector<vector<double> > mat_calor(int n, double l, double k);

vector<double> generar_vector_eccalor(int n, double l, double k, double q, double phi0);
*/

#endif