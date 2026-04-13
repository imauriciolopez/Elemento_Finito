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

template <typename T> vector<vector<T> > eye(T data, int n){
    vector<vector<T> > I(n, vector<T>(n, (T)0.0));
    for(int i=0;i<n;i++) I[i][i]=data;
    return I;
}

template <typename T> vector<vector<T> > eye(vector<T> data){
    vector<vector<T> > I(data.size(), vector<T>(data.size(), (T)0.0));
    for(int i=0;i<data.size();i++) I[i][i]=data[i];
    return I;
}

tuple<double, vector<vector<double> > > invertir_matriz(vector<vector<double> > mat, double toler=1e-12, bool det=false);

struct nodo{
    int id;
    vector<double> posicion;
    tuple<bool, double> cond_dirichlet;

    nodo(int id_=0, vector<double> pos=vector<double>{}, bool estado=false, double temp=0.0): id(id_), posicion(pos){
        cond_dirichlet=make_tuple(estado, temp);
    }

    ~nodo(){}

    void ver(){
        cout<<"id: "<<id<<", posicion: ";
        for(int j=0;j<posicion.size();j++){
            cout<<posicion[j]<<", ";
        }
        cout<<"DIRICHLET: "<<get<0>(cond_dirichlet)<<", "<<get<1>(cond_dirichlet)<<endl;
    }
};

struct material{
    double cond_neuman;
    string nombre;
    double k;

    vector<vector<double> > D;

    void ver(){
        cout<<"neuman: "<<cond_neuman<<", nombre: "<<nombre<<", k: "<<k<<endl;
        cout<<"MATIZ D:"<<endl;
        for(auto a:D){
            for(auto b:a){
                cout<<b<<", ";
            }
            cout<<endl;
        }
    }
};

class elemento{
public:
    string geometria;
    int n_dims;
    int n_nodos;
    int n_ptos_gauss_fuerza;
    int n_ptos_gauss_rigidez;
    vector<struct nodo*> nodos;
    struct material material;
    double Q;
    vector<vector<double> > nodos_posiciones;

    elemento(int dims, int nodos, int n_ptos_gauss_rigidez_, int n_ptos_gauss_fuerza_, string geo):n_dims(dims), n_nodos(nodos), geometria(geo), n_ptos_gauss_rigidez(n_ptos_gauss_rigidez_), n_ptos_gauss_fuerza(n_ptos_gauss_fuerza_){}

    virtual vector<double> Ns(const vector<double> rho) const=0;
    virtual vector<vector<double> > DNs(const vector<double> rho, bool transposicion) const=0;

    //cambiar rho para que reciba varios puntos

    //recibe un apuntador de double 
    tuple<double, vector<vector<double> >, vector<vector<double> > > crear_J(vector<double> rho) const;

    //┌virtual hace que la declaración del método exista en la clase padre, pero se tiene que definir en cada clase hijo
    //|                                                       ┌const significa que no va a cambiar el objeto original
    vector<vector<double> > matriz_rigidez_elemental() const;
    vector<double> vector_fuerza_elemental() const;

    virtual ~elemento() {}

    void ver();
};

tuple<vector<vector<double> >, vector<double> > puntos_pesos_gauss(int n, string geo);



#endif