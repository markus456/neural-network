#include "genetics.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <sstream>
#include "glew.h"
#include <SFML\Graphics.hpp>
#define GENERATIONS 1000
#define CREATURE_CAP 8
#define FOOD_CAP 8
#define SCREEN_W 1200
#define SCREEN_H 800
#define AREA_W 1000
#define AREA_H 1000
bool update();
void render();
int stringPopulation();
std::vector<std::shared_ptr<Creature>> creatures;
VectorPopulation* vpop;
sf::Time prev_click;
std::vector<sf::Vector2f> food;
sf::Vector2f mousepos;
sf::RectangleShape bg(sf::Vector2f(AREA_W,AREA_H));
sf::Shape* fToken;
sf::RenderWindow* window;
sf::Clock mouseclock,frameclock;
sf::Text* text;
sf::View camera;
int fcnt,evolution_rate,increase_mutations, decrease_mutations;
double speed, pavg, avg, current_mutation_rate,bestavg;
int counter, gen, decline_counter;
std::vector<Genome> best_so_far;
void updateCreatures();
void fastForward(int fgens)
{
	int fcounter = 0, gcounter = gen+fgens;
	std::stringstream ss;
	ss.str("");
	window->clear();
	ss << "Evolving " << fgens << " generations...";
	while(gen<gcounter) {
		avg = 0;
		for(auto& creature: creatures) {
			creature->foodPool(food);
			creature->update();
			for(auto fd = food.begin(); fd<food.end(); fd++) {
				if(distance(creature->getPos(),*fd)<15) {
					fd = food.erase(fd);
					creature->increaseFitness(1.f);
				}
			}

		}

		while(food.size()<FOOD_CAP) {
			food.push_back(sf::Vector2f(20 + rand() % (AREA_W-40), 10 + rand() % (AREA_H-20)));
		}

		if(--fcounter<0) {
			gen++;
			updateCreatures();
			fcounter = 60*evolution_rate;
		}
	}
}
void updateCreatures()
{
	srand(time(0));
	vpop->clear();
	for(int i = 0; i<creatures.size(); i++) {
		vpop->add(creatures[i]->getDNA());
	}
	vpop->update();
	auto vec = vpop->getPopulation();
	std::sort(creatures.begin(),creatures.end(),[](const std::shared_ptr<Creature>& a,const std::shared_ptr<Creature>& b) {
		return a->getFitness()>b->getFitness();
	});
	for(int i = 0; i<creatures.size(); i++) {
		creatures[i]->setDNA(vec[i]);
	}
}
int main(int argc, char** argv)
{
	fcnt = 0;
	speed = 0.5;
	gen = 1;
	evolution_rate = 10;
	counter = 300;
	bestavg = 0;
	decline_counter = 0;
	srand(time(0));
	camera.setSize(800,600);
	camera.setCenter(AREA_W/2,AREA_H/2);
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window = new sf::RenderWindow(sf::VideoMode(SCREEN_W,SCREEN_H),"Evolution 0.3",sf::Style::Default,settings);
	window->setFramerateLimit(60);
	window->setView(camera);
	bg.setOutlineColor(sf::Color::Cyan);
	bg.setOutlineThickness(4.f);
	bg.setFillColor(sf::Color(51,153,153));
	sf::CircleShape* tmp;
	tmp = new sf::CircleShape(15.f,5);
	tmp->setOutlineColor(sf::Color::Red);
	tmp->setFillColor(sf::Color(75,75,75));
	tmp->setOutlineThickness(2.f);
	tmp->setOrigin(tmp->getLocalBounds().width/2,tmp->getLocalBounds().height/2);
	fToken = tmp;
	sf::Font font;
	font.loadFromFile("revalia.ttf");
	text = new sf::Text();
	text->setFont(font);
	text->setCharacterSize(15.f);
	text->setColor(sf::Color::White);
	vpop = new VectorPopulation();
	vpop->setMutationScale(current_mutation_rate);
	for(int i = 0; i<CREATURE_CAP; i++) {
		creatures.push_back(std::shared_ptr<Creature>(new Creature(4,2,1,6)));
		creatures.back()->setPos(sf::Vector2f(10 + rand() % (AREA_W-20), 10 + rand() % (AREA_H-20)));
		creatures.back()->areaSize(AREA_W,AREA_H);
		creatures.back()->setSpeed(3);
	}
	for(int i = 0; i<FOOD_CAP; i++) {
		food.push_back(sf::Vector2f(20 + rand() % (AREA_W-40), 10 + rand() % (AREA_H-20)));
	}

	while(update()) {
		render();
	}
	return 0;
}
bool update()
{
	sf::Event event;
	while(window->pollEvent(event)) {
		switch (event.type) {
		case sf::Event::Closed:
			return false;
		case sf::Event::KeyPressed:
			switch(event.key.code) {
			case sf::Keyboard::Up:
				camera.move(0,-1);
				break;
			case sf::Keyboard::Down:
				camera.move(0,1);
				std::cout << camera.getCenter().x;
				break;
			case sf::Keyboard::Left:
				camera.move(-1,0);
				break;
			case sf::Keyboard::Right:
				camera.move(1,0);
				break;
			case sf::Keyboard::F1:
				for(int i = 0; i<10; i++) {
					updateCreatures();
				}
				break;
			case sf::Keyboard::F2:
				for(int i = 0; i<100; i++) {
					updateCreatures();
				}
				break;
			case sf::Keyboard::F3:
				fastForward(10);
				break;
			case sf::Keyboard::F4:
				fastForward(50);
				break;
			}
			break;
		case sf::Event::MouseMoved:
			if(sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
				camera.move(mousepos-window->mapPixelToCoords(sf::Mouse::getPosition(),camera));
			}
			mousepos = window->mapPixelToCoords(sf::Mouse::getPosition(),camera);
			break;
		case sf::Event::MouseWheelMoved:
			if(event.mouseWheel.delta>0) {
				camera.zoom(1.15);
			} else if(event.mouseWheel.delta<0) {
				camera.zoom(0.85);
			}
			break;
		case sf::Event::MouseButtonPressed:
			if(mouseclock.getElapsedTime().asMilliseconds()<333) {
				auto center = window->mapPixelToCoords(sf::Vector2i(event.mouseButton.x,event.mouseButton.y),camera);
				camera.setCenter(center);
				camera.setSize(sf::Vector2f(SCREEN_W,SCREEN_H));
			}
			prev_click = mouseclock.restart();
			break;
		}
	}
	double ftest = 0;
	avg = 0;
	for(auto& creature: creatures) {
		creature->foodPool(food);
		creature->update();
		for(auto fd = food.begin(); fd<food.end(); fd++) {
			if(distance(creature->getPos(),*fd)<15) {
				fd = food.erase(fd);
				creature->increaseFitness(1.f);
			}
		}
		avg += creature->getFitness();
		if(creature->getFitness()>ftest) {
			ftest = creature->getFitness();
		}
	}
	pavg = avg;
	avg = avg/creatures.size();
	if(avg>bestavg) {
		best_so_far = vpop->getPopulationRef();
	} else if(avg < bestavg*0.85) {
		vpop->revertPopulation(best_so_far);
	}

	while(food.size()<FOOD_CAP) {
		food.push_back(sf::Vector2f(20 + rand() % (AREA_W-40), 10 + rand() % (AREA_H-20)));
	}

	if(counter>1) {
		counter--;
	} else {
		updateCreatures();
		counter = 60*evolution_rate;
		std::cout << "Generation: " << gen++;
		std::cout << " fittest: " << ftest;
		std::cout << " average: " << avg << std::endl;
	}
	return true;
}
void render()
{
	std::stringstream ss;
	std::deque<Point> trails, vertices;
	window->clear();
	window->setView(camera);
	window->draw(bg);
	sf::Color color(sf::Color::Red);
	sf::ConvexShape tail;
	tail.setFillColor(sf::Color(128,0,0,128));
	for(auto& creature:creatures) {
		window->draw(*creature);
	}
	for(auto& f:food) {
		fToken->setPosition(f);
		window->draw(*fToken);
	}
	ss.str("");
	ss << "Generation: " << gen;
	text->setPosition(10,10);
	text->setString(ss.str());
	window->setView(window->getDefaultView());
	window->draw(*text);
	if(frameclock.getElapsedTime().asMilliseconds()>0) {
		ss.str("");
		ss << "FPS: " << (int)1000.f/frameclock.restart().asMilliseconds();
		text->setString(ss.str());
		text->move(0,text->getLocalBounds().height);
		window->draw(*text);
	}

	window->display();
}

int stringPopulation()
{
	Population pop;
	pop.setTarget(55.5);
	bool found = false;
	double prev = -1;
	int i = 0, gen = 0;;
	while(!found) {
		pop.update();
		if(pop.fittestGenome().second>0.999) {
			found = true;
			std::cout << "Solution found.\n";
			continue;
		}
		if(pop.fittestGenome().second>prev||i % 30000 == 0) {
			prev = pop.fittestGenome().second;
			std::cout << "\nGeneration: " << ++gen << " Target: " << pop.getTarget() << std::endl;
			std::string answer;
			for(auto& a:pop.getPopulation()) {
				answer = a.get();
				while(answer.find_first_of(':')!=std::string::npos)answer.replace(answer.find_first_of(':'),1,1,'+');
				while(answer.find_first_of(';')!=std::string::npos)answer.replace(answer.find_first_of(';'),1,1,'-');
				while(answer.find_first_of('<')!=std::string::npos)answer.replace(answer.find_first_of('<'),1,1,'*');
				while(answer.find_first_of('=')!=std::string::npos)answer.replace(answer.find_first_of('='),1,1,'/');
				std::cout << answer << std::endl;
			}
			answer = pop.fittestGenome().first;
			while(answer.find_first_of(':')!=std::string::npos)answer.replace(answer.find_first_of(':'),1,1,'+');
			while(answer.find_first_of(';')!=std::string::npos)answer.replace(answer.find_first_of(';'),1,1,'-');
			while(answer.find_first_of('<')!=std::string::npos)answer.replace(answer.find_first_of('<'),1,1,'*');
			while(answer.find_first_of('=')!=std::string::npos)answer.replace(answer.find_first_of('='),1,1,'/');
			std::cout << answer << std::endl;
			std::cout << "\nFittest individual: " << answer << " " << pop.fittestGenome().second << std::endl;
			std::cout << "Population size: " << pop.getPopulation().size() << std::endl;
			std::cin.get();
			i = 0;
		}
		gen++;

		i++;
		if(i % 1000 == 0)std::cout << ".";
	}
	std::cout << "\nGeneration: " << ++gen << " Target: " << pop.getTarget() << std::endl;
	std::string answer;
	for(auto& a:pop.getPopulation()) {
		answer = a.get();
		while(answer.find_first_of(':')!=std::string::npos)answer.replace(answer.find_first_of(':'),1,1,'+');
		while(answer.find_first_of(';')!=std::string::npos)answer.replace(answer.find_first_of(';'),1,1,'-');
		while(answer.find_first_of('<')!=std::string::npos)answer.replace(answer.find_first_of('<'),1,1,'*');
		while(answer.find_first_of('=')!=std::string::npos)answer.replace(answer.find_first_of('='),1,1,'/');
		std::cout << answer << std::endl;
	}
	answer = pop.fittestGenome().first;
	while(answer.find_first_of(':')!=std::string::npos)answer.replace(answer.find_first_of(':'),1,1,'+');
	while(answer.find_first_of(';')!=std::string::npos)answer.replace(answer.find_first_of(';'),1,1,'-');
	while(answer.find_first_of('<')!=std::string::npos)answer.replace(answer.find_first_of('<'),1,1,'*');
	while(answer.find_first_of('=')!=std::string::npos)answer.replace(answer.find_first_of('='),1,1,'/');
	std::cout << answer << std::endl;
	std::cout << "\nFittest individual: " << answer << " " << pop.fittestGenome().second << std::endl;
	std::cout << "Population size: " << pop.getPopulation().size() << std::endl;
	std::cin.get();
	return 0;
}



