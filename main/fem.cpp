#include "../include/funciones_1d.hpp"
#include "../include/utils.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <tuple>

using namespace std;

int main(){
    vector<vector<double> > D{{1.0}};
    vector<int> indices{0, 1};
    vector<vector<double> > posiciones{{0.0},
                                       {1.0}};

    int n_elems=20;
    vector<vector<double> > mat_k_global(n_elems+1, vector<double>(n_elems+1, 0.0));
    vector<double> vect_f_global(n_elems+1, 0.0);

    double l=1.0, k=10.0, cond_neumann=25.0, phi_0=5.0;

    for(int i=0;i<n_elems;i++){
        //cambiamos manualmente las condiciones para cada elemento
        //eso se haría en automático leyendo los datos directamente de GiD
        for(int j=0;j<2;j++){
            indices[j]+=1;
        }
        for(int j=0;j<2;j++){
            posiciones[j][0]+=1.0;
        }

        elem_1d_2n elemento(indices, posiciones, D, cond_neumann);

        //obtenemos la matriz de rigidez
        vector<vector<double> > mat_k=elemento.matriz_rigidez_elemental();
        //obtenemos el vector de fuerza
        vector<double> vect_f=elemento.vector_fuerza_elemental();

        
        cout<<"----------"<<(k/(l/n_elems))<<endl;
        //ensamblamos
        for(int x=0;x<2;x++){
            for(int y=0;y<2;y++){
                mat_k_global[i+x][i+y]+=mat_k[x][y]*(k/(l/n_elems));
            }
        }

        for(int j=0;j<2;j++){
            vect_f_global[i+j]+=vect_f[j];
        }
    }

    vect_f_global[0]=(cond_neumann*(l/n_elems))+(phi_0*k)/(l/n_elems);
    vect_f_global[vect_f_global.size()-1]=cond_neumann*(l/n_elems)/2.0;
    
    
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
