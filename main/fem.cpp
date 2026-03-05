#include "../include/funciones_1d.hpp"
#include "../include/utils.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <tuple>

using namespace std;

int main(){
    //OBTENCION DE DATOS
    vector<vector<double> > D{{1.0}};
    vector<int> indices{0, 1};

    int n_elems=20;
    double cond_neumann=25.0;
    double l=1.0, k=10.0, phi_0=5.0, posicion=0.0;

    //CREACION DE LOS NODOS
    vector<nodo*> nodos(n_elems+1, nullptr);
    for(int i=0;i<n_elems+1;i++){
        nodos[i]=new struct nodo(i, vector<double>{posicion}, false);
        posicion+=1.0;
    }

    get<0>(nodos[0]->cond_dirichlet)=true;
    get<1>(nodos[0]->cond_dirichlet)=15.0;
    
    get<0>(nodos[n_elems]->cond_dirichlet)=true;
    get<1>(nodos[n_elems]->cond_dirichlet)=17.0;

    //CREACION DE LOS ELEMENTOS
    vector<elem_1d_2n*> elementos(n_elems, nullptr);
    
    for(int i=0;i<n_elems;i++){
        //cambiamos manualmente las condiciones para cada elemento
        //eso se haría en automático leyendo los datos directamente de GiD
        elementos[i]=new elem_1d_2n(indices, vector<struct nodo*>{nodos[0], nodos[1]}, D, cond_neumann);
    }

    
    //ENSAMBLE
    vector<vector<double> > mat_k_global(n_elems+1, vector<double>(n_elems+1, 0.0));
    vector<double> vect_f_global(n_elems+1, 0.0);

    for(int i=0;i<n_elems;i++){
        //obtenemos la matriz de rigidez
        vector<vector<double> > mat_k=elementos[i]->matriz_rigidez_elemental();
        //obtenemos el vector de fuerza
        vector<double> vect_f=elementos[i]->vector_fuerza_elemental();

        //ensamblamos
        for(int x=0;x<2;x++) for(int y=0;y<2;y++) mat_k_global[i+x][i+y]+=mat_k[x][y];

        for(int j=0;j<2;j++) vect_f_global[i+j]+=vect_f[j];
    }

    //VERIFICACION
    for(auto a:mat_k_global){
        for(auto b:a){
            cout<<b<<", ";
        }
        cout<<endl;
    }
    cout<<"----------------"<<endl;
    for(auto a:vect_f_global){
        cout<<a<<", ";
    }
    cout<<"\n--------------"<<endl;
    cout<<"----------------"<<endl;

    //APLICACION DE CONDICIONES
    for(int i=0;i<n_elems+1;i++){
        //iterar sobre todos los nodos con temperatura fija
        if(get<0>(nodos[i]->cond_dirichlet)){
            //ajustar el vector de fuerza
            for(int j=0;j<n_elems+1;j++){
                vect_f_global[j]-=mat_k_global[j][i]*get<1>(nodos[i]->cond_dirichlet);
            }

            //anular fila y columna
            for(int j=0;j<n_elems+1;j++){
                mat_k_global[j][i]=0.0;
                mat_k_global[i][j]=0.0;
            }

            //poner 1 en la diagonal
            mat_k_global[i][i]=1.0;

            //ajustar el vector de fuerza
            vect_f_global[i]=get<1>(nodos[i]->cond_dirichlet);
        }
    }
    
    //VERIFICACION
    for(auto a:mat_k_global){
        for(auto b:a){
            cout<<b<<", ";
        }
        cout<<endl;
    }
    cout<<"----------------"<<endl;
    for(auto a:vect_f_global){
        cout<<a<<", ";
    }

    //RESOLUCION
    vector<double> phi=factorizacion_choleski_llt(mat_k_global, vect_f_global);

    //ESCRITURA
    FILE *file_temps=fopen("../archivo_temps.txt", "w");
    for(int j=0;j<phi.size()-1;j++){
        fprintf(file_temps, "%lf, ", phi[j]);
    }
    fprintf(file_temps, "%lf\n", phi[phi.size()-1]);
    fclose(file_temps);

    //LIBERACION
    for(auto e:elementos) delete e;
    for(auto n:nodos) delete n;
        
    return 0;
}

/*
double l=1.0, k=10.0, q=25.0, phi_0=5.0;

    FILE *file=fopen("archivo.txt", "w");
    FILE *file_temps=fopen("archivo_temps.txt", "w");

    for(int i=1;i<550;i++){
        vector<vector<double> > mat=mat_calor(2*i+1, l, k);
        vector<double> vect=generar_vector_eccalor(2*i+1, l, k, q, phi_0);
        vector<double> res=factorizacion_choleski_llt(mat, vect);
        fprintf(file, "%lf, %lf\n", res[ceil((2*i+1)/2)], res[res.size()-1]);

        for(int j=0;j<res.size()-1;j++){
            fprintf(file_temps, "%lf, ", res[j]);
        }
        fprintf(file_temps, "%lf\n", res[res.size()-1]);

        if(i%20==0) cout<<"N nodos: "<<i*2+1<<", avance: "<<(i*100.0)/550.0<<"%"<<endl;
    }

    fclose(file);
    fclose(file_temps);
*/
