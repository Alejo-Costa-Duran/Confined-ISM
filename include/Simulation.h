#pragma once
#include <iostream>
#ifndef SIMULATION
#define SIMULATION
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <vector>
#include <map>

class Simulation {

public:
    std::string simulation_inputs_file;
    std::map<std::string, double> parameters;
    std::string sim_type;

    Simulation();
    Simulation(std::string sim_file);
    ~Simulation(){}

    void set_par(std::string par, double value);
    void simAnnulus(std::map<std::string,double> par, std::string folder, std::string run_id, std::string fData);
    void simMosquito(std::map<std::string,double> par, std::string folder, std::string run_id, std::string fData);
    void readSimFile(std::string sim_file);
    void run(std::string sim_type, std::string folder, std::string run_id, std::string infoFile);

    void print_usage();
};

#endif