#ifndef FUNCIONES_HPP
#define FUNCIONES_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>

using namespace std;

template <typename T> vector<vector<T> > multiplicar(const vector<vector<T>>& A, const vector<vector<T>>& B){
    int n=A.size(), m=A[0].size(), p=B[0].size();

    //Verificación básica
    if(m!=B.size()){
        cout<<"DIMENSIONES INCOMPATIBLES"<<endl;
        return {};
    }
        

    vector<vector<T>> C(n, vector<T>(p, 0.0));

    for(int i=0;i<n;i++){
        for (int j=0;j<p;j++){
            for (int k=0;k<m;k++){
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    }

    return C;
}

struct nodos{
    int n_dims;
    int n_nodos;

    vector<int> indices;
    vector<vector<double> > posiciones;

    double cond_dirichlet;

    //nada conocido
    nodos(int n_dims_=0):n_nodos(0), indices(vector<int>(0)), posiciones(vector<vector<double> >(n_nodos, vector<double>(n_dims, 0.0))), n_dims(n_dims_){}

    //solo alguno de los 3 conocidos
    nodos(int n_dims_, int c_nodos):n_nodos(c_nodos), indices(vector<int>(c_nodos)), posiciones(vector<vector<double> >(c_nodos, vector<double>(3, 0.0))), n_dims(n_dims_){}
    nodos(int n_dims_, vector<int> c_indices):n_nodos(indices.size()), indices(c_indices), posiciones(vector<vector<double> >(c_indices.size(), vector<double>(3, 0.0))), n_dims(n_dims_){}
    nodos(int n_dims_, vector<vector<double> > c_posiciones):n_nodos(c_posiciones.size()), indices(vector<int>(c_posiciones.size())), posiciones(c_posiciones), n_dims(n_dims_){}

    //solo 1 desconocido
    nodos(int n_dims_, int c_nodos, vector<int> c_indices):n_nodos(c_nodos), indices(vector<int>(c_indices)), posiciones(vector<vector<double> >(c_nodos, vector<double>(3, 0.0))), n_dims(n_dims_){}
    nodos(int n_dims_, int c_nodos, vector<vector<double> > c_posiciones):n_nodos(c_nodos), indices(vector<int>(c_nodos)), posiciones(c_posiciones), n_dims(n_dims_){}
    nodos(int n_dims_, vector<int> c_indices, vector<vector<double> > c_posiciones):n_nodos(c_indices.size()), indices(c_indices), posiciones(c_posiciones), n_dims(n_dims_){}

    //ningún desconocido
    nodos(int n_dims_, int c_nodos, vector<int> c_indices, vector<vector<double> > c_posiciones):n_nodos(c_nodos), indices(c_indices), posiciones(c_posiciones), n_dims(n_dims_){}
};

struct material{
    double cond_neuman;
    string nombre;
    double k;
};

class elemento{
protected:
    int n_dims;
    int n_nodos;
    struct nodos nodos;
    struct material material;

public:
    elemento(int dims, int nodos):n_dims(dims), n_nodos(nodos) {}

    //┌virtual hace que la declaración del método exista en la clase padre, pero se tiene que definir en cada clase hijo
    //|                                                       ┌const significa que no va a cambiar el objeto original
    virtual vector<vector<double> > matriz_rigidez_elemental(const vector<double>& rho) const=0;
    virtual vector<double> vector_fuerza_elemental(const vector<double>& rho) const=0;

    virtual vector<vector<double> > Ns(const vector<double>& rho, int orden) const=0;
    virtual vector<vector<double> > DNs(const vector<double>& rho, int orden) const=0;

    /*
    vector<vector<double> > J(const vector<double>& rho){
        auto dN=DNs(rho, 1);

        int dim_nat=n_dims;  //filas
        int n_nodos=n_nodos; //columnas
        int dim=n_dims;      //dimensión física

        vector<vector<double>> J(dim, vector<double>(dim_nat, 0.0));

        for(int i=0;i<dim;i++){
            for(int j=0;j<dim_nat;j++){
                for(int a=0;a<n_nodos;a++){
                    J[i][j]+=nodos.posiciones[a][i]*dN[j][a];
                }
            }
        }

        return J;
    };*/

    virtual ~elemento() {}
};

class elem_1d_2n:public elemento{
public:
    //                    ┌1 dim
    //                    |  ┌2 nodos
    elem_1d_2n(double cond_neumann_=0.0, vector<int>indices, vector<vector<double> > posiciones):elemento(1, 2){
        nodos.n_dims=1;
        nodos.n_nodos=2;
        nodos.indices=indices;
        nodos.posiciones=posiciones;
        material.cond_neuman=cond_neumann_;
    }

    vector<vector<double> > Ns(const vector<double>& rho, int orden=1) const override{
        if(orden==1) return {{(1.0-rho[0])/2.0, (1.0+rho[0])/2.0}};
        else return {{(1.0-rho[0])/2.0}, {(1.0+rho[0])/2.0}};
    }

    vector<vector<double> > DNs(const vector<double>& rho, int orden=1) const override{
        if(orden==1) return {{-0.5, 0.5}};
        else return {{-0.5}, {0.5}};
    }

    vector<vector<double> > J(const vector<double>& rho) const{
        return {{(nodos.posiciones[1][0]-nodos.posiciones[0][0])/2.0}};
    }

    vector<vector<double> > matriz_rigidez_elemental(const vector<double>& rho) const override{
        //k = \sum_p B^T D B det|J| w_p

        //                 ┌pto de gauss
        vector<double> rho_{0};
        //       ┌peso de gauss
        double w=2.0;
        vector<vector<double> > jac=J(rho_);

        vector<vector<double> > B_T=DNs(rho_, -1);
        for(auto& a:B_T){
            for(auto& b:a){
                b/=jac[0][0];
            }
        }
        vector<vector<double> > B  =DNs(rho_);
        for(auto& a:B){
            for(auto& b:a){
                b/=jac[0][0];
            }
        }
        
        
        vector<vector<double> > k=multiplicar<double>(B_T, B);

        for(auto& a:k){
            for(auto& b:a){
                b*=w*jac[0][0]*material.k;
            }
        }

        return k;
    }

    vector<double> vector_fuerza_elemental(const vector<double>& rho) const override {
        //f = \sum_p N^T Q det|J| w_p
        
        //Punto Gauss
        vector<double> rho_={-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
        double w=1.0;

        //Funciones de forma en punto Gauss
        auto N=Ns(rho_);

        //Jacobiano
        auto jac = J(rho_);
        double Jval = jac[0][0];

        vector<double> fe(2, 0.0);

        //ciclo sobre los puntos de gauss
        for(int i=0;i<2;)
        for(int a = 0; a < 2;++a)
            fe[a]+=N[0][a]*material.cond_neuman*Jval*w;

        return fe;
    }
};

vector<double> factorizacion_choleski_llt(vector<vector<double> > arri, vector<double> vect);

vector<vector<double> > mat_calor(int n, double l, double k);

vector<double> generar_vector_eccalor(int n, double l, double k, double q, double phi0);

#endif