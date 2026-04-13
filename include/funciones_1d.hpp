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
    //                                                                                                                         ┌1 dim
    //                                                                                                                         |  ┌2 nodos
    //                                                                                                                         |  |  ┌2 ptos de integracion de funciones de forma
    //                                                                                                                         |  |  |  ┌1 pto de integracion de derivadas de funciones de forma
    elem_1d_2n(vector<int> indices, vector<struct nodo*> nodos_, vector<vector<double> > D, double cond_neumann_=0.0):elemento(1, 2, 2, 1, "lineal"){
        nodos=nodos_;
        nodos_posiciones=vector<vector<double> >(2, vector<double>(1, 0.0));
        for(int i=0;i<n_nodos;i++){
            nodos_posiciones[i]=nodos[i]->posicion;
        }
        material.D=D;
        material.cond_neuman=cond_neumann_;
    }

    vector<double> Ns(const vector<double> rho) const override;

    vector<vector<double> > DNs(const vector<double> rho, bool transposicion=true) const override;
};

class elem_1d_3n:public elemento{
public:
    //                                                                                                                         ┌1 dim
    //                                                                                                                         |  ┌3 nodos
    //                                                                                                                         |  |  ┌3 ptos de integracion de funciones de forma
    //                                                                                                                         |  |  |  ┌2 pto de integracion de derivadas de funciones de forma
    elem_1d_3n(vector<int> indices, vector<struct nodo*> nodos_, vector<vector<double> > D, double cond_neumann_=0.0):elemento(1, 3, 3, 2, "lineal"){
        nodos=nodos_;
        nodos_posiciones=vector<vector<double> >(n_nodos, vector<double>(n_dims, 0.0));
        for(int i=0;i<n_nodos;i++){
            nodos_posiciones[i]=nodos[i]->posicion;
        }
        material.D=D;
        material.cond_neuman=cond_neumann_;
    }

    vector<double> Ns(const vector<double> rho) const override;

    vector<vector<double> > DNs(const vector<double> rho, bool transposicion=true) const override;
};

class elem_2d_3n:public elemento{
public:
    //                                                                                                                         ┌2 dims
    //                                                                                                                         |  ┌3 nodos
    //                                                                                                                         |  |  ┌1 ptos de integracion de funciones de forma
    //                                                                                                                         |  |  |  ┌1 pto de integracion de derivadas de funciones de forma
    elem_2d_3n(vector<int> indices, vector<struct nodo*> nodos_, vector<vector<double> > D, double cond_neumann_=0.0):elemento(2, 3, 1, 1, "triangular"){
        nodos=nodos_;
        nodos_posiciones=vector<vector<double> >(n_nodos, vector<double>(n_dims, 0.0));
        for(int i=0;i<n_nodos;i++){
            nodos_posiciones[i]=nodos[i]->posicion;
        }
        material.D=D;
        material.cond_neuman=cond_neumann_;
    }

    vector<double> Ns(const vector<double> rho) const override;

    vector<vector<double> > DNs(const vector<double> rho, bool transposicion=true) const override;
};

class elem_3d_4n:public elemento{
public:
    //                                                                                                                         ┌3 dims
    //                                                                                                                         |  ┌4 nodos
    //                                                                                                                         |  |  ┌1 ptos de integracion de funciones de forma
    //                                                                                                                         |  |  |  ┌1 pto de integracion de derivadas de funciones de forma
    elem_3d_4n(vector<int> indices, vector<struct nodo*> nodos_, vector<vector<double> > D, double cond_neumann_=0.0):elemento(3, 4, 1, 1, "tetrahedrica"){
        nodos=nodos_;
        nodos_posiciones=vector<vector<double> >(n_nodos, vector<double>(n_dims, 0.0));
        for(int i=0;i<n_nodos;i++){
            nodos_posiciones[i]=nodos[i]->posicion;
        }
        material.D=D;
        material.cond_neuman=cond_neumann_;
    }

    vector<double> Ns(const vector<double> rho) const override;

    vector<vector<double> > DNs(const vector<double> rho, bool transposicion=true) const override;
};


/*
LEGACY
vector<vector<double> > mat_calor(int n, double l, double k);

vector<double> generar_vector_eccalor(int n, double l, double k, double q, double phi0);
*/

#endif