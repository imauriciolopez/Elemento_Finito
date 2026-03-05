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

vector<vector<double> > elem_1d_2n::matriz_rigidez_elemental() const{
    //k = \sum_p B^T D B det|J| w_p
    //                          ┌recibimos todos los puntos de gauss y sus pesos
    tuple<vector<vector<double> >, vector<double> > g_q=puntos_pesos_gauss(n_nodos-1, geometria);
    //                      ┌inicializamos la matriz k en ceros
    vector<vector<double> > k_local(n_nodos, vector<double>(n_nodos, 0.0));

    //iteramos sobre todos los puntos de gauss
    for(int punto_gauss=0;punto_gauss<n_nodos-1;punto_gauss++){
        
        //                                                                     ┌calculamos "de golpe" el jacobiano, su inversa y su determinante
        tuple<double, vector<vector<double> >, vector<vector<double> > > inf_J=crear_J(get<0>(g_q)[punto_gauss]);
        
        //                        ┌calculamos la matriz B
        //                        |                   ┌inversa del jacobiano
        //                        |                   |              ┌matriz de derivadas
        vector<vector<double> > B=multiplicar<double>(get<2>(inf_J), DNs(get<0>(g_q)[punto_gauss]));
        
        //                      ┌calculamos el producto D B
        vector<vector<double> > D_B=multiplicar<double>(material.D, B);
        //                      ┌calculamos el producto B_T D B
        vector<vector<double> > B_T_D_B=multiplicar<double>(transpuesta<double>(B), D_B);
        //                      ┌calculamos el producto B_T D B det|J| w_p
        vector<vector<double> > B_T_D_B_detJ_W=multiplicar<double>(B_T_D_B, get<0>(inf_J)*get<1>(g_q)[punto_gauss]);
        
        //sumamos punto por punto la matriz K obtenida a la matriz de rigidez elemental
        for(int i=0;i<n_nodos;i++){
            for(int j=0;j<n_nodos;j++){
                k_local[i][j]+=B_T_D_B_detJ_W[i][j];
            }
        }
        
    }
    
    //regresamos la matriz de rigidez elemental ya formada
    return k_local;
}

vector<double> elem_1d_2n::vector_fuerza_elemental() const{
    //f = \sum_p N^T Q det|J| w_p

    //                          ┌recibimos todos los puntos de gauss y sus pesos
    tuple<vector<vector<double> >, vector<double> > g_q=puntos_pesos_gauss(n_nodos, geometria);
    //                      ┌inicializamos el vector de fuerza elemental en ceros
    vector<double> f_local(n_nodos, 0.0);
    
   for(int punto_gauss=0;punto_gauss<n_nodos-1;punto_gauss++){
        
        vector<double> mat_derivadas=Ns(get<0>(g_q)[punto_gauss]);
        tuple<double, vector<vector<double> >, vector<vector<double> > > inf_J=crear_J(get<0>(g_q)[punto_gauss]);

        vector<double> N_T_Q=multiplicar<double>(mat_derivadas, material.cond_neuman*get<0>(inf_J)*get<1>(g_q)[punto_gauss]);

        //sumamos elemento por elemento la matriz K obtenida a la matriz de rigidez elemental
        for(int i=0;i<n_nodos;i++){
            f_local[i]+=N_T_Q[i];
        }
    }

    return f_local;
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