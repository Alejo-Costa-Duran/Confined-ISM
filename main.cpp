#include "manager.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>


void simCircles(int totSteps,int numRadi, double boxSize,double separation, int numAngles,double inert,
double coupling, double maxs, double timestep, double eta, double temp, std::string filename)
{
	ofstream file("./Data/Circles/"+filename);
	int seed = rand();
	Flock fl = Flock(seed,boxSize,timestep, maxs, inert, eta, temp, coupling, 1.0,0.0);
	int l = 0;
	for (int i = 1; i < numRadi + 1; i++)
		for (int j = 0; j < numAngles * i; j++)
		{
			double theta = 2.0 * M_PI * j / (1.0*numAngles * i);
			double radius = separation * i + 1.0;
			Pvector pos(radius * cos(theta), radius * sin(theta));
			Pvector vel(-maxs * sin(theta), maxs * cos(theta));
			Pvector acc(-maxs * maxs * cos(theta) / radius, -maxs * maxs * sin(theta) / radius);
			if(i==numRadi){fl.bandada.push_back(Bird(false,l, pos, vel, acc, inert, maxs,boxSize));}
			else{fl.bandada.push_back(Bird(false,l, pos, vel, acc, inert, maxs,boxSize));}
			++l;
		}
	for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, coupling, vec);
	}
	int c = 100;
	if (file.is_open())
	{
		std::string header = "Spin,SpinSq,Ener,EnerSq,TangSpd,TangSpdSq,VecX,VecXSq,VecY,VecySq,"+std::to_string(seed);
		file << header << std::endl;
		for (int i = 0; i < totSteps; i++)
		{
			fl.updateFlock(timestep);
			fl.boundary();
			if (c == 100)
			{
				cout.precision(17);
				std::vector<double> medidas = fl.measurements();
				for (unsigned int l = 0; l < medidas.size() - 1;l++) { file << std::setprecision(10) << medidas[l] << ","; }
				file << medidas[medidas.size()-1] << std::endl;
				c = 0;
				std::cout << i << "\n";
			}
			c++;
		}
		file.close();
		std::string c = "./Data/Circles/Configs/rng"+std::to_string(temp);
		char *cc =  const_cast<char*>(c.c_str());
		fl.saveFile(cc,"./Data/Circles/Configs/Enviroment"+std::to_string(temp)+
		".csv","./Data/Circles/Configs/Birds"+std::to_string(temp)+".csv",timestep);
	}
}

void simExtField(int seed,int sampling,int totSteps,int numBirds,double field, double boxSize,double inert,double coupling, double maxs, double timestep, double eta, double temp, std::string filename)
{
    ofstream filePos("./Data/Mosquito/Pos"+filename);
    ofstream fileVel("./Data/Mosquito/Vel"+filename);
    ofstream fileSpin("./Data/Mosquito/Spin"+filename);
	Flock fl = Flock(seed,boxSize,timestep, maxs, inert, eta, temp, coupling, 1.0,field);
	for(int idx=0; idx<numBirds; idx++)
	{
        double posx = 1;
		double posy = 0;
		double direction = M_PI*0.5;
		double spin = 0;
		Pvector pos(posx,posy);
		Pvector vel(0, maxs);
		Pvector acc(-spin*vel.y/inert,spin*vel.x/inert);
		fl.bandada.push_back(Bird(false,idx, pos, vel, acc, inert, maxs,boxSize));
	}
	for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, coupling, vec);
	}
	int c= sampling;
	if (filePos.is_open() && fileVel.is_open() && fileSpin.is_open())
	{
		std::string headerPos = "";
		std::string headerVel = "";
		std::string headerSpin = "";
		for(int idx=0; idx<numBirds-1;idx++)
		{
            headerPos += std::to_string(idx+1) + "PosX,"+std::to_string(idx+1)+"PosY,";
            headerVel += std::to_string(idx+1) + "VelX,"+std::to_string(idx+1)+"VelY,";
            headerSpin += std::to_string(idx+1) + "Spin,"+std::to_string(seed);
		}
		headerPos += std::to_string(numBirds)+"PosX,"+std::to_string(numBirds)+"PosY";
		headerVel += std::to_string(numBirds)+"VelX,"+std::to_string(numBirds)+"VelY";
		headerSpin += std::to_string(numBirds)+"Spin,"+std::to_string(seed);

		filePos << headerPos << std::endl;
		fileVel << headerVel << std::endl;
		fileSpin << headerSpin <<std::endl;
		for (int i = 0; i < totSteps; i++)
		{
			if (c == sampling)
			{
				cout.precision(17);
                filePos << std::setprecision(10) <<fl.bandada[0].position.x<<","<<fl.bandada[0].position.y<<"\n";
                fileVel << std::setprecision(10) <<fl.bandada[0].velocity.x<<","<<fl.bandada[0].velocity.y<<"\n";
                fileSpin <<std::setprecision(10) <<fl.bandada[0].getSpin(inert,maxs)<<"\n";
				std::cout << i << "\n";
				c=0;
			}
            fl.updateFlock(timestep);
			fl.boundary();
			c++;
		}
		filePos.close();
		fileVel.close();
		fileSpin.close();
		std::string title = "./Data/Mosquito/Configs/rng"+filename;
		char *cc =  const_cast<char*>(title.c_str());
		fl.saveFile(cc,"./Data/Mosquito/Configs/Enviroment"+filename,"./Data/Mosquito/Configs/Birds"+filename,timestep);
	}
}


void simLoops(int seed, int sampling, int totSteps, double boxSize, double inert, double timestep, double temperature, std::string filename)
{

    ofstream filePos("../Data/Mosquito/Pos"+filename);
    ofstream fileVel("../Data/Mosquito/Vel"+filename);
    ofstream fileSpin("../Data/Mosquito/Spin"+filename);

    Flock fl(seed,boxSize,timestep,1.0,inert,0.0,temperature,0.0,1.0,0.0);

    Pvector Pos(3,0);
    Pvector Vel(0,1);
    Pvector acc(-1.0/3.0,0);

    std::string headerPos = "PosX,PosY";
    std::string headerVel = "VelX,VelY";
    std::string headerSpin = "Spin";
    filePos<< headerPos <<"\n";
    fileVel << headerVel <<"\n";
    fileSpin << headerSpin <<"\n";

    fl.bandada.push_back(Bird(false,0,Pos,Vel,acc,inert,1.0,boxSize));
    int samp = sampling;
    for(int step = 0; step<totSteps;step++)
    {
        if(samp == sampling)
        {
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[0].position.x << "," <<fl.bandada[0].position.y <<"\n";
            fileVel << std::setprecision(10) << fl.bandada[0].velocity.x << "," <<fl.bandada[0].velocity.y <<"\n";
            fileSpin << std::setprecision(10) << fl.bandada[0].getSpin(inert,1.0) <<"\n";
            samp = 0;
        }
        fl.updateFlock(timestep);
        samp++;
    }
    filePos.close();
    fileVel.close();
    fileSpin.close();
}

int main()
{

    //Manager man;
    //man.Run(0.1);
    simLoops(40,100,1000000,20,1,1.0/10000.0,0,"prueba.csv");
    /*
    
    gsl_rng *rn = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rn,41);
    double boxSize=1000.0, temperature=10.0, timestep=0.001;
    int sampleRate =100;
    int  totSteps = 2000000;
    ofstream fileData("../Data/Mosquito/Info.txt");
    fileData << "boxSize,temperature,timestep,sampleRate,totSteps,inertia,kappa,spin\n";
    for(int realization=0; realization<20;realization++ )
    {
        std:vector<Flock> fl;
        double inertia = (realization+1)*50;
        std::string num_file = std::to_string(realization);
        num_file.insert(0, 2-num_file.length(), '0');
        ofstream filePos("../Data/Mosquito/Pos"+num_file+".csv");
        ofstream fileVel("../Data/Mosquito/Vel"+num_file+".csv");
        ofstream fileSpin("../Data/Mosquito/Spin"+num_file+".csv");
        std::string posHeader ="";
        std::string velHeader ="";
        std::string spinHeader="";
        for(int l=0; l<99; l++)
        {
            std::string lab = std::to_string(l);
            posHeader += lab+"PosX,"+lab+"PosY,";
            velHeader += lab+"VelX,"+lab+"VelY,";
            spinHeader += lab+"Spin,";
        }
        posHeader += "99PosX,99PosY\n";
        velHeader += "99VelX,VelPosY\n";
        spinHeader += "99Spin \n";
        filePos << posHeader;
        fileVel << velHeader;
        fileSpin << spinHeader;

        
        for(int l=0; l<100; l++)
        {
            int seed = gsl_rng_get(rn);
            fl.push_back(Flock(seed,boxSize,timestep,1.0,inertia,1.0,temperature,0.0,1.0,0.0));
            double theta=gsl_rng_uniform(rn)*2*M_PI;
            Pvector pos(2*gsl_rng_uniform(rn)-1,2*gsl_rng_uniform(rn)-1);
            Pvector vel(cos(theta),sin(theta));
            fl[l].bandada.push_back(Bird(false,0,Pvector(1,0),vel,Pvector(0,0),inertia,1.0,boxSize));
        }
        fileData<<boxSize<<","<<temperature<<","<<timestep<<","<<sampleRate<<","<<totSteps<<","<<inertia<<","<<1<<","<<0<<"\n";
        std::cout<<num_file<<"\n";

        for(int step=0; step<totSteps; step++)
        {
            if(step%sampleRate == 0)
            {
                for(int p=0; p<99; p++)
                {
                    cout.precision(17);
                    filePos<< std::setprecision(10)<< fl[p].bandada[0].position.x <<"," << fl[p].bandada[0].position.y<<",";
                    fileVel<< std::setprecision(10)<< fl[p].bandada[0].velocity.x <<"," << fl[p].bandada[0].velocity.y<<",";
                    fileSpin<< std::setprecision(10)<< fl[p].bandada[0].getSpin(inertia,1.0) <<",";
                }
                filePos<< std::setprecision(10)<< fl[99].bandada[0].position.x <<"," << fl[99].bandada[0].position.y<<"\n";
                fileVel<< std::setprecision(10)<< fl[99].bandada[0].velocity.x <<"," << fl[99].bandada[0].velocity.y<<"\n";
                fileSpin<< std::setprecision(10)<< fl[99].bandada[0].getSpin(inertia,1.0)<<"\n";
            }
            for(int l=0; l<100;l++)
            {
                fl[l].updateFlock(timestep);
                fl[l].boundary();
            }
        }
        filePos.close();
        fileVel.close();
        fileSpin.close();
    }
    fileData.close();
 */   
}
