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

// > "$1.txt" 2>&1

//CAMBIAR MATERIAL DE ELEMENTO, PASAR DE OBJETO A APUNTADOR


int main(int argc, char *argv[]){
    tuple<vector<struct nodo*>, vector<elemento*> >lista=lector_dat((string)(argv[1]));

    vector<struct nodo*> nodos=get<0>(lista);
    vector<elemento*> elementos=get<1>(lista);

    cout<<"NODOS:"<<endl;

    for(auto a:nodos){
        a->ver();
    }

    cout<<"----------------------------\nELEMENTOS:"<<endl;

    for(auto a:elementos){
        a->ver();
    }

    cout<<"----------------------------\n"<<endl;

    
    if(!nodos.empty()&&!elementos.empty()){
        cout<<"GENERACIÓN DE NODOS Y ELEMENTOS, PROCEDIENDO AL CALCULO"<<endl;
    }
    else{
        cout<<"GENERACIÓN DE NODOS Y ELEMENTOS INCORRECTA"<<endl;
        return 0;
    }

    cout<<"NODOS: "<<endl;
    for(int i=0;i<nodos.size();i++){
        nodos[i]->ver();
        cout<<"------------------------------------"<<endl;
    }
    cout<<"ELEMENTOS: "<<endl;
    for(int i=0;i<elementos.size();i++){
        elementos[i]->ver();
        cout<<"------------------------------------"<<endl;
    }

    
    //ENSAMBLE
    int n_nodos=nodos.size();
    int n_elems=elementos.size();
    vector<vector<double> > mat_k_global(n_nodos, vector<double>(n_nodos, 0.0));
    vector<double> vect_f_global(n_nodos, 0.0);


    
    cout<<"CREANDO ENSAMBLE DE MATRIZ DE RIGIDEZ Y VECTOR DE FUERZA"<<endl;
    for(int i=0;i<n_elems;i++){
        cout<<"mat_k ";
        //obtenemos la matriz de rigidez
        vector<vector<double> > mat_k=elementos[i]->matriz_rigidez_elemental();

        cout<<" vect_f ";
        //obtenemos el vector de fuerza
        vector<double> vect_f=elementos[i]->vector_fuerza_elemental();
        cout<<" finalizados "<<i<<endl;

        //ensamblamos
        for(int x=0;x<elementos[0]->n_nodos;x++){
            for(int y=0;y<elementos[0]->n_nodos;y++){
                mat_k_global[elementos[i]->nodos[x]->id][elementos[i]->nodos[y]->id]+=mat_k[x][y];
            }
            vect_f_global[elementos[i]->nodos[x]->id]+=vect_f[x];
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
        }
    }
    
    cout<<"CONDICIONES APLICADAS EXITOSAMENTE, SOLUCIONANDO SISTEMA "<<endl;

    //RESOLUCION
    vector<double> phi=factorizacion_choleski_llt(mat_k_global, vect_f_global);
    
    cout<<"SISTEMA SOLUCIONADO EXITOSAMENTE, ESCRIBIENDO RESPUESTA"<<endl;

    escribir_respuesta_estatica(phi, (string)(argv[2]));

    if(!nodos.empty()) for(int i=0;i<nodos.size();i++) delete nodos[i];
    if(!elementos.empty()) for(int i=0;i<elementos.size();i++) delete elementos[i];
    return 0;


    /*
    
     



    vector<vector<double> > D{{1.0}};
    vector<int> indices{0, 1};
    double cond_neumann=0.05;
    double l=1.0, k=10.0, phi_0=5.0, posicion=0.0, phi_1=5.0;

    FILE *file=fopen("../archivo.txt", "w");
    FILE *file_temps=fopen("../archivo_temps.txt", "w");

    for(int n_elems=1;n_elems<550;n_elems++){
        int n_nodos=(n_elems-1)*2+3;

        //CREACION DE LOS NODOS
        vector<nodo*> nodos(n_nodos, nullptr);
        for(int i=0;i<n_nodos;i++){
            nodos[i]=new struct nodo(i, vector<double>{posicion+i*(l/n_elems)}, false);
        }

        get<0>(nodos[0]->cond_dirichlet)=true;
        get<1>(nodos[0]->cond_dirichlet)=phi_0;

        //CREACION DE LOS ELEMENTOS
        vector<elem_1d_3n*> elementos(n_elems, nullptr);
        int nodo=0;
        for(int i=0;i<n_elems;i++){
            //cambiamos manualmente las condiciones para cada elemento
            //eso se haría en automático leyendo los datos directamente de GiD
            elementos[i]=new elem_1d_3n(indices, vector<struct nodo*>{nodos[nodo], nodos[nodo+1], nodos[nodo+2]}, D, cond_neumann);
            nodo+=2;
        }

        //ENSAMBLE
        vector<vector<double> > mat_k_global(n_nodos, vector<double>(n_nodos, 0.0));
        vector<double> vect_f_global(n_nodos, 0.0);

        for(int i=0;i<n_elems;i++){
            //obtenemos la matriz de rigidez
            vector<vector<double> > mat_k=elementos[i]->matriz_rigidez_elemental();
            //obtenemos el vector de fuerza
            vector<double> vect_f=elementos[i]->vector_fuerza_elemental();

            //ensamblamos
            for(int x=0;x<elementos[0]->n_nodos;x++){
                for(int y=0;y<elementos[0]->n_nodos;y++){
                    mat_k_global[elementos[i]->nodos[0]->id+x][elementos[i]->nodos[0]->id+y]+=mat_k[x][y];
                }
                vect_f_global[elementos[i]->nodos[0]->id+x]+=vect_f[x];
            }
        }

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
            }
        }

        //RESOLUCION
        vector<double> phi=factorizacion_choleski_llt(mat_k_global, vect_f_global);

        //ESCRITURA
        fprintf(file, "%lf, %lf\n", phi[ceil((n_elems)/2)], phi[phi.size()-1]);
        for(int j=0;j<phi.size()-1;j++){
            fprintf(file_temps, "%lf, ", phi[j]);
        }
        fprintf(file_temps, "%lf\n", phi[phi.size()-1]);

        //LIBERACION
        for(auto e:elementos) delete e;
        for(auto n:nodos) delete n;

        if(n_elems%20==0) cout<<"N nodos: "<<n_nodos<<", avance: "<<(n_elems*100.0)/550.0<<"%"<<endl;
    }


    fclose(file_temps);
    fclose(file);
    */
}


/*
int main(){
    //OBTENCION DE DATOS
    vector<vector<double> > D{{1.0}};
    vector<int> indices{0, 1};
    double cond_neumann=0.05;
    double l=1.0, k=10.0, phi_0=5.0, posicion=0.0, phi_1=5.0;

    FILE *file=fopen("../archivo.txt", "w");
    FILE *file_temps=fopen("../archivo_temps.txt", "w");
    
    for(int n_elems=1;n_elems<550;n_elems+=2){
        //CREACION DE LOS NODOS
        vector<nodo*> nodos(n_elems+1, nullptr);
        for(int i=0;i<n_elems+1;i++){
            nodos[i]=new struct nodo(i, vector<double>{posicion+i*(l/n_elems)}, false);
        }

        get<0>(nodos[0]->cond_dirichlet)=true;
        get<1>(nodos[0]->cond_dirichlet)=phi_0;

        //CREACION DE LOS ELEMENTOS
        vector<elem_1d_2n*> elementos(n_elems, nullptr);
        for(int i=0;i<n_elems;i++){
            //cambiamos manualmente las condiciones para cada elemento
            //eso se haría en automático leyendo los datos directamente de GiD
            elementos[i]=new elem_1d_2n(indices, vector<struct nodo*>{nodos[i], nodos[i+1]}, D, cond_neumann);
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

        //RESOLUCION
        vector<double> phi=factorizacion_choleski_llt(mat_k_global, vect_f_global);

        //ESCRITURA
        fprintf(file, "%lf, %lf\n", phi[ceil((n_elems)/2)], phi[phi.size()-1]);
        for(int j=0;j<phi.size()-1;j++){
            fprintf(file_temps, "%lf, ", phi[j]);
        }
        fprintf(file_temps, "%lf\n", phi[phi.size()-1]);
        

        //LIBERACION
        for(auto e:elementos) delete e;
        for(auto n:nodos) delete n;

        if(n_elems%20==0) cout<<"N nodos: "<<n_elems*2+1<<", avance: "<<(n_elems*100.0)/550.0<<"%"<<endl;
    }
    fclose(file_temps);
    fclose(file);
        
    return 0;
}
*/

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
