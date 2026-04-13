#ifndef GID_HPP
#define GID_HPP

#include "../include/utils.hpp"
#include "../include/funciones_1d.hpp"
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include <unordered_map>
#include <typeinfo>
#include <tuple>

using namespace std;

tuple<vector<struct nodo*>, vector<elemento*> > lector_dat(const string& dat_file, bool verbose=false);

void escribir_respuesta_estatica(vector<double> phi, string nombre);

#endif