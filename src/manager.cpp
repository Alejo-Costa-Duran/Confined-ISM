#include "Flock.h"
#include "manager.h"
#include "SFML/Window.hpp"
#include "SFML/Graphics.hpp"
#include <math.h>

#define BOID_AMOUNT 1
#define Pi 3.14159265359


Manager::Manager()
{
	//this->boidsSize = rand() % 10 - 3;
	sf::VideoMode desktop = sf::VideoMode::getDesktopMode();
	this->_window_width = 800;
	this->_window_height = 800;
	this->_window.create(sf::VideoMode(sf::Vector2u(_window_width, _window_height), desktop.bitsPerPixel), "Boids", sf::Style::None);
	printInstructions();
}


void Manager::Run(float fric)
{
	double time = 1.0 / 1000;
	srand(5);
	double maxs=1.0;
	double Temp = 0.0;
	double inert = 1;
	double coupling = 0;
	double friction = 0;
	flock = Flock(0,10.0,time, maxs, inert, friction, Temp,coupling, 1.0,1.0);
    int c=0;
	for (int i = 1; i < 1+1; i++)
    {
		for (int j = 0; j < 1*i; j++)
		{
            double theta = 2.0 * M_PI * (j / (10.0 * i));
			double radius = 5*i;//radius = 0.6 * i + 2.0;
			//double posX = 0.9*i;
			//double posY = 0.9*j;
			Pvector pos(radius*cos(theta),radius*sin(theta));
			Pvector vel(-sin(theta),cos(theta));
			Pvector acc(-maxs * maxs * cos(theta) / radius, -maxs * maxs * sin(theta) / radius);
			createBoid(false, c, pos, vel, acc, inert, maxs, sf::Color::Green, sf::Color::Green);
			c += 1;
			
		}
		}

    /*
    for(int idx=0; idx<100; idx++)
	{
        double posx = 2*(gsl_rng_uniform(flock.m_mt)-0.5);
		double posy = 2*(gsl_rng_uniform(flock.m_mt)-0.5);
		double direction = 0*2*Pi*gsl_rng_uniform(flock.m_mt);
		double spin = 0*(gsl_rng_uniform(flock.m_mt)-0.5);
		Pvector pos(posx,posy);
		Pvector vel(maxs * cos(direction), maxs * sin(direction));
		Pvector acc(-spin*vel.y/inert,spin*vel.x/inert);
		createBoid(false,idx, pos, vel, acc, inert, maxs, sf::Color::Green, sf::Color::Green);
	}
	for(int l=101;l<111;l++)
    {createBoid(true,l,Pvector(5+2.0*l/100.0,0),Pvector(0,1),Pvector(0,0),inert,maxs,sf::Color::Green, sf::Color::Green);}
	*/
	for (unsigned int idx = 0; idx < flock.bandada.size(); idx++)
	{
		std::vector<Bird> vec = flock.flocking(idx);
		flock.bandada[idx].calcForce(inert, coupling, vec,0);
	}

	//Whole block of text can probably simplified in a function as well in order to remove redundancy
	sf::Font font;

	sf::Text fpsText("Frames per Second: ", font);
	fpsText.setFillColor(sf::Color::White);
	fpsText.setCharacterSize(12);
	fpsText.setPosition(sf::Vector2f(_window_width - 162, 0));
	_window.setFramerateLimit(6000);
	// Clock added to calculate frame rate, may cause a small amount of slowdown?
	sf::Clock clock;
	//clock.restart();
	std::cout << flock.c1 << "\n" << flock.c0 << "\n" << flock.c2 << "\n";
	std::cout<< flock.bandada.size()<<"\n";
	int i = 0;
	while (_window.isOpen()) {
		clock.restart();
		float fps = 1 / time; // 1 / refresh time = estimate of fps
		//for(auto& b : flock.bandada){std::cout<<b.lambda<<"\t";}
		HandleInput();
		Render(fpsText, fps, 1.0 / 1000);
		time = clock.getElapsedTime().asSeconds();
		//std::cout<< flock.totalSpin()<<"\t"<<i<<"\n";
		++i;
		std::cout<<flock.bandada[0].position.getNorm()<<"\n";
		//std::cout << flock.bandada[0].velocity.dotProd(flock.bandada[0].acc) << "\n";
		//std::cout << i << "\n";
		//++i;
	}
}

//Method of passing text needs refactoring
void Manager::Render(sf::Text fpsText, double fps, double time)
{
	double lengthX = 30.0;
	double lengthY = 30.0;
	_window.clear();
	//fpsText.setString("Frames per Second: " + to_string(int(1.0/time + 0.5)));
	_window.draw(fpsText);
	sf::CircleShape circ;
	circ.setPointCount(50);
	circ.setPosition({_window_width/2-80,_window_height/2-80});
	circ.setRadius(80);
	circ.setFillColor(sf::Color(150,50,250,256));
	circ.setOutlineColor(sf::Color(150, 50, 250));
	circ.setOutlineThickness(1);

	sf::CircleShape circ1;
	circ1.setPointCount(50);
	circ1.setPosition({_window_width/2-360,_window_height/2-360});
	circ1.setRadius(360);
	circ1.setFillColor(sf::Color(150,50,250,256));
	circ1.setOutlineColor(sf::Color(150, 50, 250));
	circ1.setOutlineThickness(1);


	float r = _window_height/20.0;
	Pvector centerMass(0,0);
	for(Bird b : flock.bandada){centerMass.addVector(b.position);}
	centerMass.x = centerMass.x/1100;
	centerMass.y = centerMass.y/1100;
	_window.draw(circ);
	_window.draw(circ1);

	// Draws all of the Boids out, and applies functions that are needed to update.

	flock.updateFlock(time,0);
	//flock.boundary();
	for (unsigned int i = 0; i < shapes.size(); i++) {
		//Pvector vel = flock.bandada[i].velocity;
		//Pvector pos = flock.bandada[i].position;
		//float x = vel.radius*cos(pos.theta)-vel.theta*sin(pos.theta);
		//float y = vel.radius*sin(pos.theta)+vel.theta*cos(pos.theta);
		Pvector pos = flock.bandada[i].polarToScreen(_window_width, _window_height, flock.boxSize, flock.boxSize);
		shapes[i].setPosition(sf::Vector2f(pos.x,pos.y));
		shapes[i].setOrigin(sf::Vector2f(1, 1));
		_window.draw(shapes[i]);
	}
	/*
	circ.setPosition(shapes[0].getPosition().x-r,shapes[0].getPosition().y-r);
	_window.draw(circ);
	for(int l=0; l<lengthX*2;l++)
	{
	for(int j = 0; j<lengthY*2;j++)
	{
		float y = l*_window_height/20.0;
		float x1 = j*_window_width/20.0;
		sf::Vertex line[] =	{sf::Vertex(sf::Vector2f(x1, y)),sf::Vertex(sf::Vector2f(x1+100.0, y))};
		sf::Vertex line2[] =	{sf::Vertex(sf::Vector2f(y, x1)),sf::Vertex(sf::Vector2f(y, x1+100.0))};
	_window.draw(line, 2, sf::Lines);
	_window.draw(line2, 2, sf::Lines);
	}
	}
	*/
	_window.display();
	//std::cout<<flock.flocking(1099).size()<<"\t"<<flock.bandada[1099].boxX<<"\n";

}

void Manager::HandleInput()
{

}


void Manager::createBoid(bool fix,int idx, Pvector pos, Pvector vel, Pvector acc, double inertia, double maxs, sf::Color fillColor, sf::Color outlineColor)
{
	int size = 1;
	Bird b(fix,idx, pos, vel, acc, inertia, maxs,30.0);
	float lengthX = 30.0, lengthY = 30.0;
	sf::CircleShape shape(size, 3);
	shape.setPosition(sf::Vector2f(b.polarToScreen(_window_width, _window_height, lengthX, lengthY).x, b.polarToScreen(_window_width, _window_height, lengthX, lengthY).y));
	shape.setFillColor(fillColor);
	shape.setOutlineColor(outlineColor);
	shape.setOutlineThickness(.5);

	flock.bandada.push_back(b);
	shapes.push_back(shape);

	// New Shape is drawn
	_window.draw(shapes[shapes.size() - 1]);
}

void Manager::printInstructions()
{
	std::cout << string(100, '\n');

}
