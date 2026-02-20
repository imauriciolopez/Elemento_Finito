#include "funciones_fem.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main(){
    vector<vector<double> > v_1{{1.0, 2.0, 3.0}};
    vector<vector<double> > v_2{{4.0}};

    vector<vector<double> > v_3=multiplicar<double>(v_2, v_1);

    for(auto a:v_3){
        for(auto b:a){
            cout<<b<<", ";
        }
        cout<<endl;
    }
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
