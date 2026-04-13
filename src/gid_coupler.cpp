#include "../include/utils.hpp"
#include "../include/funciones_1d.hpp"
#include "../include/gid_coupler.hpp"

#include <tuple>
#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <omp.h>

tuple<vector<struct nodo*>, vector<elemento*> > lector_dat(const string& dat_file, bool verbose){
    ifstream archivo(dat_file.c_str());

    if(!archivo){
        cout<<"NO SE ENCONTRÓ EL ARCHIVO"<<endl;
        return make_tuple(vector<struct nodo*>{}, vector<elemento*>{});
    }

    int dimensiones, n_nodos_p_elemento, n_elems, n_mats, n_nodos, n_conds_dirichlet, n_conds_neuman_lin, n_conds_neuman_plan, n_conds_neuman_vol, buff_i;
    double buff_d;
    string line, buff_string;
    bool inverse;

    getline(archivo, line);                //# Datos generales ################################
    getline(archivo, line);                //# Problem_name
    getline(archivo, line);                //Calor
    getline(archivo, line);                //# Solver_type
    getline(archivo, line);                //Cholesky
    getline(archivo, line);                //# Tolerance
    getline(archivo, line);                //1e-5
    getline(archivo, line);                //# Max_steps
    getline(archivo, line);                //10000
    getline(archivo, line);                //# Gauss_points_quantity
    getline(archivo, line);                //1 (esto está hardcodeado en cada tipo de elemento)
    getline(archivo, line);                //# Inverse_strategy
    getline(archivo, line);                //(estrategia de inversa para 2D shells)
    istringstream k(line);                 
    string t_i;
    k>>t_i;
    if(t_i=="Pseudo_inverse") inverse=false;
    else inverse=true;
    getline(archivo, line);                //# Malla ################################
    getline(archivo, line);                //# Dimension
    getline(archivo, line);                //(n dimensiones)
    istringstream l(line);                 
    l>>dimensiones;
    if(dimensiones<1||dimensiones>3){
        cout<<"LA GEOMETRÍA NO ES COMPATIBLE CON EL PROBLEMA "<<dimensiones<<endl;
        return make_tuple(vector<struct nodo*>{}, vector<elemento*>{});
    }
    getline(archivo, line);                //# Nodos por elemento
    getline(archivo, line);                //(n nodos por elemento)
    istringstream l_1(line);                 
    l_1>>n_nodos_p_elemento;
    getline(archivo, line);                //# Numero de elementos
    getline(archivo, line);                //(n elementos)
    istringstream l_2(line);                 
    l_2>>n_elems;
    getline(archivo, line);                //# Materiales ###########################
    getline(archivo, line);                //# Numero de materiales
    getline(archivo, line);                //(n materiales)
    istringstream l_3(line);                 
    l_3>>n_mats;
    vector<vector<double> > materiales(n_mats, vector<double>(4, 0.0));
    getline(archivo, line);                //# Listado de materiales
    for(int i=0;i<n_mats;i++){
        getline(archivo, line);            //(n_mat, conductividad térmica, densidad, calor específico, M)
        istringstream l_4(line);
        l_4>>buff_i;
        for(int j=0;j<4;j++){
            l_4>>buff_d;
            materiales[i][j]=buff_d;
        }
    }
    getline(archivo, line);                //# Nodos ################################
    getline(archivo, line);                //# Numero de nodos
    getline(archivo, line);                //(n nodos)
    istringstream l_5(line);                 
    l_5>>n_nodos;
    getline(archivo, line);                //# Listado de nodos
    vector<vector<double> > posiciones_nodos(n_nodos, vector<double>(dimensiones, 0.0));
    for(int i=0;i<n_nodos;i++){
        getline(archivo, line);            //(n_nodo, X, Y, Z)
        istringstream l_6(line);
        l_6>>buff_i;
        for(int j=0;j<dimensiones;j++){
            l_6>>buff_d;
            posiciones_nodos[i][j]=buff_d;
        }
    }
    
    getline(archivo, line);                //# Listado de elementos
    vector<vector<int> > nodos_de_elementos(n_elems, vector<int>(n_nodos_p_elemento+1, -1));
    for(int i=0;i<n_elems;i++){
        getline(archivo, line);            //(n_elem, nodo_1, nodo_2, ...)
        istringstream l_7(line);
        l_7>>buff_i;
        vector<string> buff_string_vect;
        while(l_7>>buff_string){
            buff_string_vect.push_back(buff_string);
        }

        if     (buff_string_vect[buff_string_vect.size()-1]=="Linear")        nodos_de_elementos[i][n_nodos_p_elemento]=0;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Triangle")      nodos_de_elementos[i][n_nodos_p_elemento]=1;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Quadrilateral") nodos_de_elementos[i][n_nodos_p_elemento]=2;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Circle")        nodos_de_elementos[i][n_nodos_p_elemento]=3;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Tetrahedra")    nodos_de_elementos[i][n_nodos_p_elemento]=4;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Hexahedra")     nodos_de_elementos[i][n_nodos_p_elemento]=5;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Prism")         nodos_de_elementos[i][n_nodos_p_elemento]=6;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Sphere")        nodos_de_elementos[i][n_nodos_p_elemento]=7;
        else if(buff_string_vect[buff_string_vect.size()-1]=="Points")        nodos_de_elementos[i][n_nodos_p_elemento]=8;
        else                                                                  nodos_de_elementos[i][n_nodos_p_elemento]=-1;

        buff_string_vect.pop_back();

        for(int j=0;j<buff_string_vect.size();j++){
            nodos_de_elementos[i][j]=stoi(buff_string_vect[j]);
        }

    }
    getline(archivo, line);                //# Materiales de elementos
    vector<int> materiales_por_elemento(n_elems, 0);
    for(int i=0;i<n_elems;i++){
        getline(archivo, line);            //(n_elem, material)
        istringstream l_8(line);
        l_8>>buff_i;
        l_8>>buff_i;
        materiales_por_elemento[i]=buff_i;
    }
    getline(archivo, line);                //# Temperaturas #########################
    getline(archivo, line);                //# Número de condiciones
    getline(archivo, line);                //(n conds)
    istringstream l_9(line);                 
    l_9>>n_conds_dirichlet;
    vector<tuple<int, double> > conds_dirichlet(n_conds_dirichlet);
    getline(archivo, line);                //# Listado de condiciones
    for(int i=0;i<n_conds_dirichlet;i++){
        getline(archivo, line);            //(n_nodo, material)
        istringstream l_10(line);
        l_10>>buff_i>>buff_d;
        get<0>(conds_dirichlet[i])=buff_i;
        get<1>(conds_dirichlet[i])=buff_d;
    }
    getline(archivo, line);                //# Terminos fuente ######################
    getline(archivo, line);                //# Superficie lineal
    getline(archivo, line);                //# Número de condiciones
    getline(archivo, line);                //(n conds)
    istringstream l_11(line);                 
    l_11>>n_conds_neuman_lin;
    vector<tuple<int, double> > conds_neuman_lin(n_conds_neuman_lin);
    vector<double> conds_neuman_lin_comp(n_elems, 0.0);
    getline(archivo, line);                //# Listado de condiciones
    getline(archivo, line);                //# Element Heat
    for(int i=0;i<n_conds_neuman_lin;i++){
        getline(archivo, line);            //(n_elem, calor)
        istringstream l_12(line);
        l_12>>buff_i>>buff_d;
        get<0>(conds_neuman_lin[i])=buff_i;
        get<1>(conds_neuman_lin[i])=buff_d;
        conds_neuman_lin_comp[buff_i-1]=buff_d;
    }
    getline(archivo, line);                //# Superficie planar
    getline(archivo, line);                //# Número de condiciones
    getline(archivo, line);                //(n conds)
    istringstream l_13(line);                 
    l_13>>n_conds_neuman_plan;
    vector<tuple<int, double> > conds_neuman_plan(n_conds_neuman_plan);
    vector<double> conds_neuman_plan_comp(n_elems, 0.0);
    getline(archivo, line);                //# Listado de condiciones
    getline(archivo, line);                //# Element Heat
    for(int i=0;i<n_conds_neuman_plan;i++){
        getline(archivo, line);            //(n_elem, calor)
        istringstream l_14(line);
        l_14>>buff_i>>buff_d;
        get<0>(conds_neuman_plan[i])=buff_i;
        get<1>(conds_neuman_plan[i])=buff_d;
        conds_neuman_plan_comp[buff_i-1]=buff_d;
    }

    getline(archivo, line);                //# Fuente volumétrica
    getline(archivo, line);                //# Número de condiciones
    getline(archivo, line);                //(n conds)
    istringstream l_15(line);                 
    l_15>>n_conds_neuman_vol;
    vector<tuple<int, double> > conds_neuman_vol(n_conds_neuman_vol);
    vector<double> conds_neuman_vol_comp(n_elems, 0.0);
    getline(archivo, line);                //# Listado de condiciones
    getline(archivo, line);                //# Element Heat
    for(int i=0;i<n_conds_neuman_vol;i++){
        getline(archivo, line);            //(n_elem, calor)
        istringstream l_16(line);
        l_16>>buff_i>>buff_d;
        get<0>(conds_neuman_vol[i])=buff_i;
        get<1>(conds_neuman_vol[i])=buff_d;
        conds_neuman_vol_comp[buff_i-1]=buff_d;
    }
    
    //---------------------------------------------------------------------
    //TERMINA LA LECTURA DEL ARCHIVO, EMPIEZA LA CREACIÓN DE LOS ELEMENTOS
    //---------------------------------------------------------------------

    if(verbose==true){
        cout<<"Dims: "<<dimensiones<<endl;
        cout<<"Nodos por elemento: "<<n_nodos_p_elemento<<endl;
        cout<<"Numero de elementos: "<<n_elems<<endl;
        cout<<"Numero de nodos: "<<n_nodos<<endl;
        cout<<"Numero de materiales: "<<n_mats<<endl;
        cout<<"Numero conds dirichlet: "<<n_conds_dirichlet<<endl;
        cout<<"Numero conds neumann lin: "<<n_conds_neuman_lin<<endl;
        cout<<"Numero conds neumann plan: "<<n_conds_neuman_plan<<endl;
        cout<<"Materiales: "<<endl;
        for(auto a:materiales){
            for(auto b:a){
                cout<<b<<", ";
            }
            cout<<endl;
        }
        cout<<"Nodos:"<<endl;
        for(auto a:posiciones_nodos){
            for(auto b:a){
                cout<<b<<", ";
            }
            cout<<endl;
        }
        cout<<"Elementos:"<<endl;
        for(auto a:nodos_de_elementos){
            for(auto b:a){
                cout<<b<<", ";
            }
            cout<<endl;
        }
        cout<<"Materiales por elemento:"<<endl;
        for(auto a:materiales_por_elemento){
            cout<<a<<", ";
        }
        cout<<"\nTemps nodos:"<<endl;
        for(int i=0;i<conds_dirichlet.size();i++){
            cout<<get<0>(conds_dirichlet[i])<<", "<<get<1>(conds_dirichlet[i])<<endl;
        }
        cout<<"Calor lin:"<<endl;
        for(int i=0;i<conds_neuman_lin.size();i++){
            cout<<get<0>(conds_neuman_lin[i])<<", "<<get<1>(conds_neuman_lin[i])<<endl;
        }
        cout<<"Calor planar:"<<endl;
        for(int i=0;i<conds_neuman_plan.size();i++){
            cout<<get<0>(conds_neuman_plan[i])<<", "<<get<1>(conds_neuman_plan[i])<<endl;
        }
    }
    
    
    //Creando los nodos (sin corrección para 1D -> sin proyección en barra)
    vector<nodo*> nodos(n_nodos, nullptr);
    for(int i=0;i<n_nodos;i++){
        nodos[i]=new struct nodo(i, posiciones_nodos[i]);
    }
    //aplicación de condiciones de dirichlet a nodos
    for(int i=0;i<n_conds_dirichlet;i++){
        get<0>(nodos[get<0>(conds_dirichlet[i])-1]->cond_dirichlet)=true;
        get<1>(nodos[get<0>(conds_dirichlet[i])-1]->cond_dirichlet)=get<1>(conds_dirichlet[i]);
    }


    //Creando los elementos (cualquier geometría)
    vector<elemento*> elementos(n_elems, nullptr);
    for(int i=0;i<n_elems;i++){
        //caso elemento lineal
        if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==0){
            if(nodos_de_elementos[i].size()<3){
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }

            //elementos 2 nodos
            else if(nodos_de_elementos[i].size()==3){
                nodos_de_elementos[i].pop_back();
                elementos[i]=new elem_1d_2n(nodos_de_elementos[i], 
                                            vector<struct nodo*>{nodos[nodos_de_elementos[i][0]-1], 
                                                                 nodos[nodos_de_elementos[i][1]-1]
                                            }, 
                                            eye<double>(materiales[materiales_por_elemento[i]-1][0], dimensiones), 
                                            conds_neuman_lin_comp[i]);
            }
            //elementos 3 nodos
            else if(nodos_de_elementos[i].size()==4){
                nodos_de_elementos[i].pop_back();
                elementos[i]=new elem_1d_3n(nodos_de_elementos[i], 
                                            vector<struct nodo*>{nodos[nodos_de_elementos[i][0]-1], 
                                                                 nodos[nodos_de_elementos[i][1]-1], 
                                                                 nodos[nodos_de_elementos[i][2]-1]
                                            }, 
                                            eye<double>(materiales[materiales_por_elemento[i]-1][0], dimensiones), 
                                            conds_neuman_lin_comp[i]);
            }

            else{
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }
        }
        //caso elemento triangular
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==1){
            if(nodos_de_elementos[i].size()<4){
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }
            //elementos 3 nodos
            else if(nodos_de_elementos[i].size()==4){
                nodos_de_elementos[i].pop_back();
                elementos[i]=new elem_2d_3n(nodos_de_elementos[i], 
                                            vector<struct nodo*>{nodos[nodos_de_elementos[i][0]-1], 
                                                                 nodos[nodos_de_elementos[i][1]-1],
                                                                 nodos[nodos_de_elementos[i][2]-1]
                                            }, 
                                            eye<double>(materiales[materiales_por_elemento[i]-1][0], dimensiones), 
                                            conds_neuman_plan_comp[i]);
            }
            else{
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }
        }
        //caso elemento cuadrático
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==2){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS CUADRADOS"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento circular
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==3){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS CIRCULARES"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento Tetrahedra
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==4){
            if(nodos_de_elementos[i].size()<5){
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }
            //elementos 3 nodos
            else if(nodos_de_elementos[i].size()==5){
                nodos_de_elementos[i].pop_back();
                elementos[i]=new elem_3d_4n(nodos_de_elementos[i], 
                                            vector<struct nodo*>{nodos[nodos_de_elementos[i][0]-1], 
                                                                 nodos[nodos_de_elementos[i][1]-1],
                                                                 nodos[nodos_de_elementos[i][2]-1],
                                                                 nodos[nodos_de_elementos[i][3]-1]
                                            }, 
                                            eye<double>(materiales[materiales_por_elemento[i]-1][0], dimensiones), 
                                            conds_neuman_vol_comp[i]);
            }
            else{
                cout<<"NO SE CAPTURARON LOS NODOS CORRECTAMENTE"<<endl;
                return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
            }
        }
        //caso elemento Hexahedra
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==2){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS HEXAHÉDRICOS"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento Prism
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==2){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS PRISMÁTICOS"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento Sphere
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==2){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS ESFÉRICOS"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento Points
        else if(nodos_de_elementos[i][nodos_de_elementos[i].size()-1]==2){
            cout<<"NO HAY IMPLEMENTACIÓN DE ELEMENTOS POINTS"<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
        //caso elemento no definido
        else{
            cout<<"NO SE CAPTURÓ CORRECTAMENTE LA GEOMETRÍA DEL ELEMENTO "<<i<<endl;
            return make_tuple(vector<struct nodo*>{}, vector<elemento*>{}); 
        }
    }
    
    return make_tuple(nodos, elementos);

    //return make_tuple(vector<struct nodo*>{}, vector<elemento*>{});
}

void escribir_respuesta_estatica(vector<double> phi, string nombre){
    string nombre_1=nombre+".post.res";
    ofstream file(nombre_1.c_str());

    if(!file.is_open()){
        cout << "No se pudo abrir el archivo "<<nombre_1<<endl;
        return;
    }

    file<<"GiD Post Results File 1.0"<<endl;
    file<<"Result \"Calor\" \"Heat\" 1 Scalar OnNodes"<<endl;
    file<<"Values"<<endl;

    for(int i=0;i<phi.size();i++){
        file<<"      "<<i+1<<"     "<<scientific<<phi[i]<<endl;
    }
    file<<"End Values"<<endl;

    string nombre_2=nombre+".results";
    ofstream file_2(nombre_2.c_str());

    if(!file_2.is_open()){
        cout << "No se pudo abrir el archivo "<<nombre_2<<endl;
        return;
    }

    file_2<<"GiD Post Results File 1.0"<<endl;
    file_2<<"Result \"Calor\" \"Heat\" 1 Scalar OnNodes"<<endl;
    file_2<<"Values"<<endl;

    for(int i=0;i<phi.size();i++){
        file_2<<"      "<<i+1<<"     "<<scientific<<phi[i]<<endl;
    }
    file_2<<"End Values"<<endl;
}