#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include <unordered_map>
#include <typeinfo>
#include <tuple>

using namespace std;

void hello_world_utils();

vector<double> factorizacion_choleski_llt(vector<vector<double> > arri, vector<double> vect);

template <typename T> vector<vector<T> > multiplicar(const vector<vector<T>> A, const vector<vector<T>> B){
    int n=A.size(), m=A[0].size(), p=B[0].size();

    //verificación básica
    if(m!=B.size()){
        throw runtime_error("Error: las dimensiones no coinciden.");
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

template <typename T> vector<vector<T> > multiplicar(vector<vector<T>> A, T elem){
    int n=A.size(), m=A[0].size();

    vector<vector<T>> C(n, vector<T>(m, 0.0));

    for(int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            C[i][j]+=A[i][j]*elem;
        }
    }

    return C;
}

template <typename T> vector<T> multiplicar(vector<T> A, T elem){
    int n=A.size();

    vector<T> C(n, 0.0);

    for(int i=0;i<n;i++) C[i]+=A[i]*elem;

    return C;
}

template <typename T> vector<vector<T> > transpuesta(const vector<vector<T>> A){
    int n=A.size(), m=A[0].size();

    vector<vector<T>> C(m, vector<T>(n, 0.0));

    for(int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            C[j][i]=A[i][j];
        }
    }

    return C;
}

tuple<double, vector<vector<double> > > invertir_matriz(vector<vector<double> > mat, double toler=1e-12);

struct nodo{
    int id;
    vector<double> posicion;
    tuple<bool, double> cond_dirichlet;

    nodo(int id=0, vector<double> pos=vector<double>{}, bool estado=false, double temp=0.0): posicion(pos){
        cond_dirichlet=make_tuple(estado, temp);
    }

    ~nodo(){}
};

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

    vector<vector<double> > D;
};

class elemento{
public:
    string geometria;
    int n_dims;
    int n_nodos;
    vector<struct nodo*> nodos;
    struct material material;
    double Q;
    vector<vector<double> > nodos_posiciones;


    elemento(int dims, int nodos, string geo):n_dims(dims), n_nodos(nodos), geometria(geo){}

    //┌virtual hace que la declaración del método exista en la clase padre, pero se tiene que definir en cada clase hijo
    //|                                                       ┌const significa que no va a cambiar el objeto original
    virtual vector<vector<double> > matriz_rigidez_elemental() const=0;
    virtual vector<double> vector_fuerza_elemental() const=0;

    virtual vector<double> Ns(const vector<double> rho) const=0;
    virtual vector<vector<double> > DNs(const vector<double> rho, bool transposicion) const=0;

    //cambiar rho para que reciba varios puntos

    //recibe un apuntador de double 
    tuple<double, vector<vector<double> >, vector<vector<double> > > crear_J(vector<double> rho) const{
        vector<vector<double> > J=multiplicar<double>(DNs(rho, true), nodos_posiciones);

        tuple<double, vector<vector<double> > > det_J_inv=invertir_matriz(J);

        return make_tuple(get<0>(det_J_inv), J, get<1>(det_J_inv));
    }

    virtual ~elemento() {}
};

tuple<vector<vector<double> >, vector<double> > puntos_pesos_gauss(int n, string geo);

#endif