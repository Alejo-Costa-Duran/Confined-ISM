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
	this->_window.create(sf::VideoMode(_window_width, _window_height, desktop.bitsPerPixel), "Boids", sf::Style::None);
	printInstructions();
}


void Manager::Run(float fric)
{
	double time = 1.0 / 100;
	srand(5);
	double maxs = 1;
	double inert = 5.;
	double coupling = 1.0;
	flock = Flock(0,10.0,time, maxs, inert, 0.000, 0.0, coupling, 1.0);
	int c = 0;
	for (int i = 1; i < 10 + 1; i++)
		for (int j = 0; j < i * 20; j++)
		{
			//double theta = 2*Pi* static_cast <double> (rand()) / static_cast <double> (RAND_MAX);;
			double theta = 2.0 * Pi * (j / (20.0 * i));
			double radius = 0.8 * i + 2.0;
			Pvector pos(radius * cos(theta), radius * sin(theta));
			Pvector vel(-maxs * sin(theta), maxs * cos(theta));
			Pvector acc(-maxs * maxs * cos(theta) / radius, -maxs * maxs * sin(theta) / radius);
			createBoid(c, pos, vel, acc, inert, maxs, sf::Color::Green, sf::Color::Green);
			c += 1;
		}
	for (unsigned int idx = 0; idx < flock.bandada.size(); idx++)
	{
		std::vector<Bird> vec = flock.flocking(idx);
		flock.bandada[idx].calcForce(inert, coupling, vec);
	}
	//Whole block of text can probably simplified in a function as well in order to remove redundancy
	sf::Font font;
	font.loadFromFile("consola.ttf");


	sf::Text fpsText("Frames per Second: ", font);
	fpsText.setFillColor(sf::Color::White);
	fpsText.setCharacterSize(12);
	fpsText.setPosition(_window_width - 162, 0);
	_window.setFramerateLimit(1);
	// Clock added to calculate frame rate, may cause a small amount of slowdown?
	sf::Clock clock;
	//clock.restart();
	std::cout << flock.c1 << "\n" << flock.c0 << "\n" << flock.c2 << "\n";
	int i = 0;
	while (_window.isOpen()) {
		clock.restart();
		float fps = 1 / time; // 1 / refresh time = estimate of fps
		//for(auto& b : flock.bandada){std::cout<<b.lambda<<"\t";}
		HandleInput();
		Render(fpsText, fps, 1.0 / 100);
		time = clock.getElapsedTime().asSeconds();
		//std::cout << flock.bandada[0].velocity.dotProd(flock.bandada[0].acc) << "\n";
		//std::cout << i << "\n";
		//++i;
	}
}

//Method of passing text needs refactoring
void Manager::Render(sf::Text fpsText, double fps, double time)
{
	double lengthX = 10.0;
	double lengthY = 10.0;
	_window.clear();
	//fpsText.setString("Frames per Second: " + to_string(int(1.0/time + 0.5)));
	_window.draw(fpsText);
	sf::CircleShape circ;
	float r = _window_height/20.0;
	circ.setRadius(_window_height/20.0);
	circ.setOrigin(0, 0);
	circ.setOutlineThickness(1.0);
	circ.setFillColor(sf::Color(0,0,0,0));
	circ.setOutlineColor(sf::Color::White);

	// Draws all of the Boids out, and applies functions that are needed to update.

	flock.updateFlock(time);
	flock.boundary();
	for (unsigned int i = 0; i < shapes.size(); i++) {
		//Pvector vel = flock.bandada[i].velocity;
		//Pvector pos = flock.bandada[i].position;
		//float x = vel.radius*cos(pos.theta)-vel.theta*sin(pos.theta);
		//float y = vel.radius*sin(pos.theta)+vel.theta*cos(pos.theta);
		Pvector pos = flock.bandada[i].polarToScreen(_window_width, _window_height, flock.boxSize, flock.boxSize);
		shapes[i].setPosition(pos.x, pos.y);
		shapes[i].setRotation(-atan2(flock.bandada[i].velocity.y, flock.bandada[i].velocity.x) * 180 / 3.14159 + 90);
		shapes[i].setOrigin(1, 1);
		shapes[i].setFillColor(sf::Color::Green);
		_window.draw(shapes[i]);
	}
	circ.setPosition(shapes[1099].getPosition().x-r,shapes[1099].getPosition().y-r);
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
	_window.display();
	//std::cout<<flock.flocking(1099).size()<<"\t"<<flock.bandada[1099].boxX<<"\n"; 

}

void Manager::HandleInput()
{
	sf::Event event;
	while (_window.pollEvent(event)) {
		// "close requested" event: we close the window
		// Implemented alternate ways to close the window. (Pressing the escape, X, and BackSpace key also close the program.)
		if ((event.type == sf::Event::Closed) ||
			(event.type == sf::Event::KeyPressed &&
				event.key.code == sf::Keyboard::Escape) ||
			(event.type == sf::Event::KeyPressed &&
				event.key.code == sf::Keyboard::BackSpace) ||
			(event.type == sf::Event::KeyPressed &&
				event.key.code == sf::Keyboard::X))
		{
			_window.close();
		}
		if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::P)
		{
			std::cout << flock.bandada.size() << "\n";
		}
	}
}


void Manager::createBoid(int idx, Pvector pos, Pvector vel, Pvector acc, double inertia, double maxs, sf::Color fillColor, sf::Color outlineColor)
{
	int size = 1;
	Bird b(idx, pos, vel, acc, inertia, maxs,10.0);
	float lengthX = 30, lengthY = 30;
	sf::CircleShape shape(size, 3);
	shape.setPosition(b.polarToScreen(_window_width, _window_height, lengthX, lengthY).x, b.polarToScreen(_window_width, _window_height, lengthX, lengthY).y);
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
