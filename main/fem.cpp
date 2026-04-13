#include "../include/funciones_1d.hpp"
#include "../include/utils.hpp"
#include "../include/gid_coupler.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <tuple>
#include <fstream>
#include <unistd.h>

using namespace std;

//CAMBIAR MATERIAL DE ELEMENTO, PASAR DE OBJETO A APUNTADOR

int main(int argc, char *argv[]){
    tuple<vector<struct nodo*>, vector<elemento*> >lista=lector_dat((string)(argv[1]), true);

    vector<struct nodo*> nodos=get<0>(lista);
    vector<elemento*> elementos=get<1>(lista);

    if(!nodos.empty()&&!elementos.empty()){
        cout<<"GENERACIÓN DE NODOS Y ELEMENTOS, PROCEDIENDO AL CALCULO"<<endl;
    }
    else{
        cout<<"ERROR DETECTADO ";
        if(nodos.empty()){
            cout<<"EN LA CREACIÓN DE NODOS, SALIENDO DEL PROGRAMA"<<endl;
        }
        else if(nodos.empty()&&elementos.empty()){
            cout<<"EN LA CREACIÓN DE NODOS Y ELEMENTOS, SALIENDO DEL PROGRAMA"<<endl;
        }
        else{
            cout<<"EN LA CREACIÓN DE ELEMENTOS, SALIENDO DEL PROGRAMA"<<endl;
        }
        return 0;
    }

    int n_nodos=nodos.size();
    int n_elems=elementos.size();

    if(n_nodos<100){
        cout<<"NODOS:"<<endl;
        for(auto a:nodos){
            a->ver();
        }
        cout<<"----------------------------\nELEMENTOS:"<<endl;

        for(auto a:elementos){
            a->ver();
        }
        cout<<"----------------------------\n"<<endl;
    }
    
    //ENSAMBLE
    vector<vector<double> > mat_k_global(n_nodos, vector<double>(n_nodos, 0.0));
    vector<double> vect_f_global(n_nodos, 0.0);

    cout<<"CREANDO ENSAMBLE DE MATRIZ DE RIGIDEZ Y VECTOR DE FUERZA"<<endl;

    for(int i=0;i<n_elems;i++){
        if(n_nodos<100)cout<<"mat_k ";
        
        //obtenemos la matriz de rigidez
        vector<vector<double> > mat_k=elementos[i]->matriz_rigidez_elemental();

        if(n_nodos<100)cout<<" vect_f ";
        //obtenemos el vector de fuerza
        vector<double> vect_f=elementos[i]->vector_fuerza_elemental();
        if(n_nodos<100)cout<<" finalizados "<<i<<endl;

        //ensamblamos
        for(int x=0;x<elementos[0]->n_nodos;x++){
            for(int y=0;y<elementos[0]->n_nodos;y++){
                mat_k_global[elementos[i]->nodos[x]->id][elementos[i]->nodos[y]->id]+=mat_k[x][y];
            }
            vect_f_global[elementos[i]->nodos[x]->id]+=vect_f[x];
        }
        if(n_nodos>1000){
            if(i%int(floor(n_elems/10))==0) cout<<"AVANCE: "<<(n_elems/i)*100<<"%"<<endl;
        }
    }
    
    cout<<"ENSAMBLE CREADO EXITOSAMENTE, APLICANDO CONDICIONES"<<endl;

    //APLICACION DE CONDICIONES
    for(int i=0;i<n_nodos;i++){
        //iterar sobre todos los nodos con temperatura fija
        if(get<0>(nodos[i]->cond_dirichlet)){
            //ajustar el vector de fuerza
            for(int j=0;j<n_nodos;j++){
                vect_f_global[j]-=mat_k_global[j][i]*get<1>(nodos[i]->cond_dirichlet);
            }

            //anular fila y columna
            for(int j=0;j<n_nodos;j++){
                mat_k_global[j][i]=0.0;
                mat_k_global[i][j]=0.0;
            }

            //poner 1 en la diagonal
            mat_k_global[i][i]=1.0;

            //ajustar el vector de fuerza
            vect_f_global[i]=get<1>(nodos[i]->cond_dirichlet);

            if(n_nodos<100) cout<<"NODO "<<i<<", temp: "<<get<1>(nodos[i]->cond_dirichlet)<<endl;
        }
    }

    cout<<"CONDICIONES APLICADAS EXITOSAMENTE, SOLUCIONANDO SISTEMA "<<endl;

    if(n_nodos<100){
        for(auto a:mat_k_global){
            for(auto b:a){
                cout<<b<<", ";
            }
            cout<<endl;
        }
        cout<<"----------------------"<<endl;
        for(auto a:vect_f_global){
            cout<<a<<", ";
        }
        cout<<"----------------------"<<endl;
    }
    
    //RESOLUCION
    vector<double> phi=factorizacion_choleski_llt(mat_k_global, vect_f_global);
    
    cout<<"SISTEMA SOLUCIONADO EXITOSAMENTE, ESCRIBIENDO RESPUESTA"<<endl;
    
    if(n_nodos<100){
        for(auto a:phi){
            cout<<a<<", ";
        }
        cout<<"----------------------"<<endl;
    }
    
    escribir_respuesta_estatica(phi, (string)(argv[2]));

    if(!nodos.empty()) for(int i=0;i<nodos.size();i++) delete nodos[i];
    if(!elementos.empty()) for(int i=0;i<elementos.size();i++) delete elementos[i];
    return 0;
}
