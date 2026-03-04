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
    elem_1d_2n(vector<int> indices, vector<vector<double> > posiciones, vector<vector<double> > D, double cond_neumann_=0.0):elemento(1, 2, "lineal"){
        nodos.n_dims=1;
        nodos.n_nodos=2;
        nodos.indices=indices;
        nodos.posiciones=posiciones;
        material.cond_neuman=cond_neumann_;
        material.D=D;
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