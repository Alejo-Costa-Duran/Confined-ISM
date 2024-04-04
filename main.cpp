#include "manager.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>
#include <Simulation.h>
#include <chrono>
#include <algorithm>

int main(int arg_count, char* args[])
{
   if(arg_count>1)
   {
    std::string folder = args[1];
    std::string run_id;
    std::string data_file = args[3];
    std::cout << "Los archivos se van a guardar en la carpeta: " + folder + "\n";
    std::cout << "Con el sufijo: " + run_id<<"\n";
    std::cout << "Los parámetros de la simulación se tomaron del archivo" + data_file<<"\n";
    Simulation s(data_file);
    int counter =0;
    std::vector<double> inertias ={0.0003,0.005,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,5};
    std::vector<double> fields = {0.0,0.001,0.005,0.01,0.1,0.25,0.5,0.75,1.0,1.05,1.5,1.75,2.0};
    for(auto chi: inertias)
    {
        for(auto f: fields)
        {        
            std::cout << "Amount of runs: " << inertias.size()*fields.size() <<"\n";
            int seed = rand();
            s.set_par("Inertia",chi);
            s.set_par("Field",f);
            s.set_par("Seed",seed);
            run_id = std::to_string(counter);
            unsigned int number_of_zeros = 3 - run_id.length();
            run_id.insert(0, number_of_zeros, '0');
            counter +=1;
            auto start = std::chrono::high_resolution_clock::now();
            s.simMosquito(s.parameters, folder, run_id,"/Info");
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            std::cout << "\n Finished simulation on file id: " + run_id;
            std::cout << "\n It took: " << duration.count() << " seconds \n";
        }
    }
   }
   else
   {
    std::cout << "Error: Es necesario proveer el nombre de la carpeta y el sufijo de los archivos \n";    return 1;
   }
}
