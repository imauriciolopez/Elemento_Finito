#include "../include/utils.hpp"

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
#include <omp.h>

void hello_world_utils(){
    cout<<"Hello, world! (utils)"<<endl;
}

vector<double> factorizacion_choleski_llt(vector<vector<double> > arri, vector<double> vect) {
    int i, j, k, n, n_elements=arri.size();
    double suma;

    vector<vector<double> > mat(n_elements, vector<double>(n_elements, 0.0));
    vector<double> be(n_elements, 0.0);
    vector<double> salida(n_elements, 0.0);
    
    for(i=0;i<n_elements;i++){
        for(j=0;j<=i;j++){
            suma=0.0;
            for(k=0;k<j;k++){
                suma+=mat[i][k]*mat[j][k];
            }
            if(i==j){
                if(arri[i][i]-suma<0){
                    cout<<"EEEERRRRROR"<<endl;
                }
                mat[i][j]=sqrt(arri[i][i]-suma);
            }
            else{
                mat[i][j]=(arri[i][j]-suma)/mat[j][j];
            }
        }
    }

    for(i=0;i<n_elements;i++){
        suma=0.0;
        for(j=0;j<i;j++){
            suma+=mat[i][j]*be[j];
        }
        be[i]=(vect[i]-suma)/mat[i][i];
    }

    for(i=n_elements-1;i>=0;i--){
        suma=0.0;
        for(j=i+1;j<n_elements;j++){
            suma+=mat[j][i]*salida[j];
        }
        salida[i]=(be[i]-suma)/mat[i][i];
    }

    return salida;
}

tuple<vector<vector<double> >, vector<double> > puntos_pesos_gauss(int n, string geo){
    //FORMATO:                     PUNTO               PESO
    if(geo=="lineal"){
        switch(n){
            case(1):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{0.0               }};
                vector<double>          elem_2{               2.0               };
                return make_tuple(elem_1, elem_2);
            }
            case(2):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.5773502691896257}, 
                                                vector<double>{-0.5773502691896257}};
                //elem 2
                vector<double>          elem_2{ 1.0               , 
                                                1.0               };
                return make_tuple(elem_1, elem_2);
            }
            case(3):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.0000000000000000},
                                                vector<double>{ 0.7745966692414834}, 
                                                vector<double>{-0.7745966692414834}};
                //elem 2
                vector<double> elem_2{ 0.8888888888888888,
                                        0.5555555555555556, 
                                        0.5555555555555556};
                return make_tuple(elem_1, elem_2);
            }
            case(4):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.3399810435848563}, 
                                                vector<double>{-0.3399810435848563},  
                                                vector<double>{ 0.8611363115940526}, 
                                                vector<double>{-0.8611363115940526}};
                //elem 2
                vector<double> elem_2{0.6521451548625461, 
                                        0.6521451548625461,
                                        0.3478548451374538, 
                                        0.3478548451374538};
                return make_tuple(elem_1, elem_2);
            }
            case(5):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.0               },  
                                                vector<double>{ 0.5384693101056831}, 
                                                vector<double>{-0.5384693101056831},  
                                                vector<double>{ 0.9061798459386640}, 
                                                vector<double>{-0.9061798459386640}};
                //elem 2
                vector<double> elem_2{0.5688888888888889, 
                                        0.4786286704993665,
                                        0.4786286704993665,
                                        0.2369268850561891,
                                        0.2369268850561891};
                return make_tuple(elem_1, elem_2);
            }
            case(6):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.6612093864662645},
                                                vector<double>{-0.6612093864662645},
                                                vector<double>{ 0.2386191860831969},
                                                vector<double>{-0.2386191860831969},
                                                vector<double>{ 0.9324695142031521},
                                                vector<double>{-0.9324695142031521}};
                //elem 2
                vector<double> elem_2{ 0.3607615730481386,
                                        0.3607615730481386,
                                        0.4679139345726910,
                                        0.4679139345726910,
                                        0.1713244923791704,
                                        0.1713244923791704};
                return make_tuple(elem_1, elem_2);
            }
            case(7):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.0               },
                                                vector<double>{ 0.4058451513773972},
                                                vector<double>{-0.4058451513773972},
                                                vector<double>{ 0.7415311855993945},
                                                vector<double>{-0.7415311855993945},
                                                vector<double>{ 0.9491079123427585},
                                                vector<double>{-0.9491079123427585}};
                //elem 2
                vector<double> elem_2{ 0.4179591836734694,
                                        0.3818300505051189,
                                        0.3818300505051189,
                                        0.2797053914892766,
                                        0.2797053914892766,
                                        0.1294849661688697,
                                        0.1294849661688697};
                return make_tuple(elem_1, elem_2);
            }
            case(8):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.1834346424956498},
                                                vector<double>{-0.1834346424956498},
                                                vector<double>{ 0.5255324099163290},
                                                vector<double>{-0.5255324099163290},
                                                vector<double>{ 0.7966664774136267},
                                                vector<double>{-0.7966664774136267},
                                                vector<double>{ 0.9602898564975363},
                                                vector<double>{-0.9602898564975363}};
                //elem 2
                vector<double> elem_2{ 0.3626837833783620,
                                        0.3626837833783620,
                                        0.3137066458778873,
                                        0.3137066458778873,
                                        0.2223810344533745,
                                        0.2223810344533745,
                                        0.1012285362903763,
                                        0.1012285362903763};
                return make_tuple(elem_1, elem_2);
            }
            case(9):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.0000000000000000},
                                                vector<double>{ 0.8360311073266358},
                                                vector<double>{-0.8360311073266358},
                                                vector<double>{ 0.9681602395076261},
                                                vector<double>{-0.9681602395076261},
                                                vector<double>{ 0.3242534234038089},
                                                vector<double>{-0.3242534234038089},
                                                vector<double>{ 0.6133714327005904},
                                                vector<double>{-0.6133714327005904}};
                //elem 2
                vector<double> elem_2{ 0.3302393550012598,
                                        0.1806481606948574,
                                        0.1806481606948574,
                                        0.0812743883615744,
                                        0.0812743883615744,
                                        0.3123470770400029,
                                        0.3123470770400029,
                                        0.2606106964029354,
                                        0.2606106964029354};
                return make_tuple(elem_1, elem_2);
            }
            case(10):{
                //elem 1
                vector<vector<double> > elem_1{vector<double>{ 0.1488743389816312},
                                                vector<double>{-0.1488743389816312},
                                                vector<double>{ 0.4333953941292472},
                                                vector<double>{-0.4333953941292472},
                                                vector<double>{ 0.6794095682990244},
                                                vector<double>{-0.6794095682990244},
                                                vector<double>{ 0.8650633666889845},
                                                vector<double>{-0.8650633666889845},
                                                vector<double>{ 0.9739065285171717},
                                                vector<double>{-0.9739065285171717}};
                //elem 2
                vector<double> elem_2{ 0.2955242247147529,
                                        0.2955242247147529,
                                        0.2692667193099963,
                                        0.2692667193099963,
                                        0.2190863625159820,
                                        0.2190863625159820,
                                        0.1494513491505806,
                                        0.1494513491505806,
                                        0.0666713443086881,
                                        0.0666713443086881};
                return make_tuple(elem_1, elem_2);
            }
            default:{return make_tuple(vector<vector<double> >{{}}, vector<double>{});}
        }
    }
    else if(geo=="triangular"){
        return make_tuple(vector<vector<double> >{{}}, vector<double>{});
    }
    else if(geo=="tetrahedrica"){
        return make_tuple(vector<vector<double> >{{}}, vector<double>{});
    }
    else if(geo=="piramidal"){
        return make_tuple(vector<vector<double> >{{}}, vector<double>{});
    }
    else{
        return make_tuple(vector<vector<double> >{{}}, vector<double>{});
    }

}

tuple<double, vector<vector<double> > > invertir_matriz(vector<vector<double> > mat, double toler){
    //caso matriz cuadrada
    int n=mat.size(), m=mat[0].size();
    if(n>3){
        return make_tuple(0.0, vector<vector<double> >{{}});
    }

    if(n==m){
        if(n==1){
            return make_tuple(mat[0][0], vector<vector<double> > {{1.0/mat[0][0]}});
        }
        else if(n==2){
            double det=mat[0][0]*mat[1][1]
                      -mat[0][1]*mat[1][0]
                      +toler;//para evitar divisiones entre 0

            return make_tuple(det-toler, vector<vector<double> > {{ mat[1][1]/det, -mat[0][1]/det},
                                                                  {-mat[1][0]/det,  mat[0][0]/det}});
        }
        else{
            double det=mat[0][0]*mat[1][1]*mat[2][2]
                      +mat[0][1]*mat[1][2]*mat[2][0]
                      +mat[0][2]*mat[1][0]*mat[2][1]

                      -mat[0][2]*mat[1][1]*mat[2][0]
                      -mat[0][0]*mat[1][2]*mat[2][1]
                      -mat[0][1]*mat[1][0]*mat[2][2]
                      +toler;//para evitar divisiones entre 0

            return make_tuple(det-toler, vector<vector<double> > {{(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])/det, (mat[2][1]*mat[0][2]-mat[0][1]*mat[2][2])/det, (mat[0][1]*mat[1][2]-mat[1][1]*mat[0][2])/det},
                                                                  {(mat[2][0]*mat[1][2]-mat[1][0]*mat[2][2])/det, (mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2])/det, (mat[1][0]*mat[0][2]-mat[0][0]*mat[1][2])/det},
                                                                  {(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1])/det, (mat[2][0]*mat[0][1]-mat[0][0]*mat[2][1])/det, (mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1])/det}});
        }
    }
    //caso matriz no cuadrada, calculamos pseudo-inversa de moore-penrose

    //esto es util cuando el jacobiano no es cuadrado, y esto pasa cua do se hace una proyección de 2 espacios físicos difenretes
    //por ejemplo, si adaptamos una malla 2D a la superficie de una esfera, el jacobiano sería de tamaño 2x3

    else{
        //caso mas filas que columnas
        //se calcula:
        //A^+= (A_T A)^-1 A^T
        if(n>m){
            return make_tuple(0.0, 
                                multiplicar<double>(
                                    get<1>(invertir_matriz(
                                        multiplicar<double>(
                                            transpuesta<double>(mat), 
                                            mat
                                        )
                                    )), 
                                    transpuesta<double>(mat)
                                )
                            );
        }
        //caso mas columnas que filas
        //se calcula:
        //A^+= A^T (A A_T)^-1 
        else{
            return make_tuple(0.0, 
                                multiplicar<double>(
                                    transpuesta<double>(mat), 
                                    get<1>(invertir_matriz(
                                        multiplicar<double>(mat, 
                                            transpuesta<double>(mat)
                                        )
                                    ))
                                )
                            );
        }
    }
}

