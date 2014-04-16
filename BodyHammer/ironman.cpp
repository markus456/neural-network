#include "genetics.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <sstream>
//#include "glew.h"
#include <SFML\Graphics.hpp>
#define GENERATIONS 1000
#define SCREEN_W 1200
#define SCREEN_H 800
#define X_RES 16
#define Y_RES 16
bool update();
void render();
std::vector<sf::RectangleShape> input_grid;
sf::RenderWindow* window;
sf::Text* text;
Beholder* beholder;
int fcnt,evolution_rate,increase_mutations, decrease_mutations;
double speed, pavg, avg, current_mutation_rate;
int counter, gen;

int wmain(void) {	
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window = new sf::RenderWindow(sf::VideoMode(SCREEN_W,SCREEN_H),"The Beholder",sf::Style::Default,settings);
	beholder = new Beholder(X_RES,Y_RES,1,(X_RES*Y_RES)/2);
	sf::Font font;
	font.loadFromFile("revalia.ttf");
	text = new sf::Text();
	text->setFont(font);
	text->setCharacterSize(24.f);
	text->setColor(sf::Color::White);
	for(int y = 0;y<Y_RES;y++){
		for(int x = 0;x<X_RES;x++){
			input_grid.push_back(sf::RectangleShape(sf::Vector2f(32.f,32.f)));
			input_grid.back().setPosition(32.f*x,32.f*y);
			input_grid.back().setFillColor(sf::Color::White);
		}
	}
	text->setPosition(0,Y_RES*32+32);
	while(update())render();
	return 0;
}
bool update() {
	sf::Event event;
	int ans = -1;
	std::stringstream ss;
	sf::String tmp;
	while(window->pollEvent(event)) {

		switch (event.type) {
		case sf::Event::Closed:
			return false;
		case sf::Event::KeyPressed:
			if(event.key.code==sf::Keyboard::F1){
				for(auto& a:input_grid){
					a.setFillColor(sf::Color::White);
				}
			}else if(event.key.code==sf::Keyboard::F2){
				ans = beholder->analyze(input_grid);
				std::stringstream ss;
				ss << "The Beholder sees a: " << ans;
				text->setString(sf::String(ss.str()));
			}
			break;
		case sf::Event::TextEntered:
				tmp = event.text.unicode;
				beholder->train(input_grid,tmp.toAnsiString()[0]);
				ss << "The Beholder is trained with: " << tmp.toAnsiString() << " Global error: " << beholder->getGlobalError();
				text->setString(ss.str());
			break;
		case sf::Event::MouseMoved:
			if(sf::Mouse::isButtonPressed(sf::Mouse::Left)){
			for(auto& a:input_grid){
				if(a.getLocalBounds().contains(sf::Vector2f(event.mouseMove.x,event.mouseMove.y)-a.getPosition())){
					a.setFillColor(sf::Color::Black);
					break;
				}
			}
			}
			break;
		}
	}
	return true;
}
void render() {
	window->clear();
	for(auto& a:input_grid){
		window->draw(a);
	}
	window->draw(*text);
	window->display();
}