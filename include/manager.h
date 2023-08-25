#include "Flock.h"
#include "SFML/Window.hpp"
#include "SFML/Graphics.hpp"

#ifndef MANAGER_H
#define MANAGER_H

class Manager {
private:
	sf::RenderWindow _window;
	int _window_width;
	int _window_height;

	Flock flock;
	float boidsSize;
	vector<sf::CircleShape> shapes;
	void Render(sf::Text text, double fps, double time);
	void createBoid(bool fix,int idx, Pvector pos, Pvector vel, Pvector acc, double inertia, double maxs, sf::Color fillColor, sf::Color outlineColor);
	void HandleInput();

public:
	// For console instructions
	static void printInstructions();
	Manager();
	void Run(float fric);
	void boundary(double rmin, double rmax);
};


#endif
