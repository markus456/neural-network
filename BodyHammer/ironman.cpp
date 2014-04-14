#include "genetics.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <sstream>
//#include "glew.h"
#include <SFML\Graphics.hpp>
#define GENERATIONS 1000
#define CREATURE_CAP 15
#define ZONES 15
#define SCREEN_W 1200
#define SCREEN_H 800
bool mainLoop();

std::vector<std::shared_ptr<Creature>> creatures;
VectorPopulation* vpop;
std::vector<Point> zones;
std::vector<sf::RectangleShape> zoneTokens;
sf::Shape* cToken;
sf::Shape* fToken;
sf::RenderWindow* window;
sf::Clock* frameclock;
sf::Text* text;
int fcnt,evolution_rate,increase_mutations, decrease_mutations;
double speed, pavg, avg, current_mutation_rate;
int counter, gen;
void updateCreatures();
void fastForward(int fgens) {
	text->setCharacterSize(30);	
	int fcounter = 0, gcounter = gen+fgens;	
	std::stringstream ss;
	while(gen<gcounter){
		ss.str("");
		window->clear();
		ss << "Evolving " << fgens << " generations...";
		text->setString(ss.str());
		text->setPosition(sf::Vector2f(SCREEN_W/2-text->getLocalBounds().width/2,SCREEN_H/2-text->getLocalBounds().height/2));
		window->draw(*text);
		ss.str("");
		ss << "Current generation: " << gen;
		text->setString(ss.str());
		text->setPosition(sf::Vector2f(SCREEN_W/2-text->getLocalBounds().width/2,SCREEN_H/2+text->getLocalBounds().height/2));
		window->draw(*text);
		window->display();
		for(auto& creature: creatures) {
			creature->zoneAt(zones.front());
			creature->update();
			for(auto& z:zones){
				if(z.isInside(creature->getPoint())) {
					creature->increaseFitness(1.f/60.f);
				}
			}
		}
		if(fcounter>1) {
			fcounter--;
		} else {
			gen++;
			updateCreatures();
			fcounter = 60*evolution_rate;
		}
	}
}
void updateCreatures() {
	srand(time(0));
	vpop->clear();
	for(int i = 0; i<creatures.size(); i++) {
		vpop->add(creatures[i]->getDNA());
	}
	vpop->update();
	auto vec = vpop->getPopulation();
	std::sort(creatures.begin(),creatures.end(),[](const std::shared_ptr<Creature>& a,const std::shared_ptr<Creature>& b){return a->getFitness()>b->getFitness();});
	for(int i = 0; i<creatures.size(); i++) {
		creatures[i]->setDNA(vec[i]);
	}
}
int wmain(void) {
	fcnt = 0;
	speed = 0.5;
	current_mutation_rate = 0.1;
	gen = 1;
	increase_mutations = decrease_mutations = 0;
	evolution_rate = 2;
	counter = 300;
	srand(time(0));
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window = new sf::RenderWindow(sf::VideoMode(SCREEN_W,SCREEN_H),"Evolution 0.3",sf::Style::Default,settings);
	sf::CircleShape* tmp = new sf::CircleShape(15.f,3);
	tmp->setOutlineColor(sf::Color::Cyan);
	tmp->setFillColor(sf::Color(75,75,75));
	tmp->setOutlineThickness(1.f);
	tmp->setOrigin(tmp->getLocalBounds().width/2,tmp->getLocalBounds().height/2);
	cToken = tmp;
	tmp = new sf::CircleShape(15.f,5);
	tmp->setOutlineColor(sf::Color::Red);
	tmp->setFillColor(sf::Color(75,75,75));
	tmp->setOutlineThickness(2.f);
	tmp->setOrigin(tmp->getLocalBounds().width/2,tmp->getLocalBounds().height/2);
	fToken = tmp;
	frameclock = new sf::Clock();
	sf::Font font;
	font.loadFromFile("revalia.ttf");
	text = new sf::Text();
	text->setFont(font);
	text->setCharacterSize(15.f);
	text->setColor(sf::Color::White);
	vpop = new VectorPopulation();
	vpop->setMutationScale(current_mutation_rate);
	for(int i = 0; i<CREATURE_CAP; i++) {
		creatures.push_back(std::shared_ptr<Creature>(new Creature(5,2,1,6)));
		creatures.back()->setPos(10 + rand() % (SCREEN_W-20), 10 + rand() % (SCREEN_H-20));
		creatures.back()->areaSize(SCREEN_W,SCREEN_H);
		creatures.back()->setSpeed(2);
	}
	for(int i = 0;i<ZONES;i++){
		Point p( rand() % (SCREEN_W - 150),rand() % (SCREEN_H - 150) ,100 + rand() % 50, 100 + rand() % 50);
		bool inside = true;
		while(zones.size()>0&&inside){
			inside = false;
			p = Point( rand() % (SCREEN_W - 150),rand() % (SCREEN_H - 150) ,100 + rand() % 50, 100 + rand() % 50);
			for(auto& zn:zones){
				if(zn.isInside(p)){
					inside = true;
					break;
				}
			}
		}
		zones.push_back(p);
		zoneTokens.push_back(sf::RectangleShape(sf::Vector2f(zones.back().w(),zones.back().h())));
		zoneTokens.back().setPosition(zones.back().x(),zones.back().y());
		zoneTokens.back().setOutlineThickness(2.f);
		zoneTokens.back().setOutlineColor(sf::Color::Cyan);
		zoneTokens.back().setFillColor(sf::Color(0,0,255,128));
	}

	while(mainLoop());
	return 0;
}
bool update() {
	sf::Event event;
	while(window->pollEvent(event)) {
		switch (event.type) {
		case sf::Event::Closed:
			return false;
		case sf::Event::KeyPressed:
			switch(event.key.code) {
			case sf::Keyboard::F1:
				if(evolution_rate>1) {
					evolution_rate--;
					speed += 3;
				}
				for(auto& creature: creatures) {
					creature->setSpeed(speed);
				}
				std::cout << "Generation length: " << evolution_rate << std::endl;
				break;
			case sf::Keyboard::F2:
				if(speed>3) {
					evolution_rate++;
					speed -= 3;
				}
				for(auto& creature: creatures) {
					creature->setSpeed(speed);
				}
				std::cout << "Generation length: " << evolution_rate << std::endl;
				break;
			case sf::Keyboard::F3:

				fastForward(50);
				break;
			case sf::Keyboard::F4:
				fastForward(250);
				break;
			}
			break;
		case sf::Event::MouseButtonPressed:

			break;
		}
	}
	double ftest = 0;
	avg = 0;

	for(auto& creature: creatures) {
		std::sort(zones.begin(),zones.end(),[&](const Point& a,const Point& b) {
			return a.distance(creature->getPoint())<b.distance(creature->getPoint());
		});
		creature->zoneAt(zones.front());
		creature->update();
		creature->setInside(false);
		for(auto& z:zones){
			if(z.isInside(creature->getPoint())) {
				creature->setInside(true);
				creature->increaseFitness(-creature->getFitness()*0.01);
				break;
			}
		}
		avg += creature->getFitness();
		if(creature->getFitness()>ftest){
			ftest = creature->getFitness();
		}
	}
	pavg = avg;
	avg = avg/creatures.size();

	if(counter>1) {
		counter--;
	} else {
		if(gen % 5==0){
			zones.clear();
			zoneTokens.clear();
			for(int i = 0;i<ZONES;i++){
				Point p( rand() % (SCREEN_W - 150),rand() % (SCREEN_H - 150) ,100 + rand() % 50, 100 + rand() % 50);
				bool inside = true;
				while(zones.size()>0&&inside){
					inside = false;
					p = Point( rand() % (SCREEN_W - 150),rand() % (SCREEN_H - 150) ,100 + rand() % 50, 100 + rand() % 50);
					for(auto& zn:zones){
						if(zn.isInside(p)){
							inside = true;
							break;
						}
					}
				}
				zones.push_back(p);
				zoneTokens.push_back(sf::RectangleShape(sf::Vector2f(zones.back().w(),zones.back().h())));
				zoneTokens.back().setPosition(zones.back().x(),zones.back().y());
				zoneTokens.back().setOutlineThickness(2.f);
				zoneTokens.back().setOutlineColor(sf::Color::Cyan);
				zoneTokens.back().setFillColor(sf::Color(0,0,255,128));
			}
		}
		for(auto& creature:creatures){
				creature->increaseFitness(1);

		}
		updateCreatures();
		counter = 60*evolution_rate;
		std::cout << "Generation: " << gen++;
		std::cout << " fittest: " << ftest;
		std::cout << " average: " << avg << std::endl;
	}
	return true;
}
void render() {
	std::stringstream ss;
	std::deque<Point> trails;
	window->clear();
	ss << "Generation: " << gen;
	text->setPosition(10,10);
	text->setString(ss.str());
	window->draw(*text);
	ss.str("");
	std::sort(creatures.begin(),creatures.end(),[](const std::shared_ptr<Creature>& a,const std::shared_ptr<Creature>& b){
		return a->getFitness()>b->getFitness();
	});
	sf::Color color(180,200,255);
	for(auto& creature:creatures){
		trails = creature->getTrail();
		double trailw = Creature::TRAIL_LENGTH+1;
		for(auto a = trails.begin();a!=trails.end();a++){
			if(a+1!=trails.end()){
				double ang = (*a).angle(*(a+1));
				sf::RectangleShape tail(sf::Vector2f((*a).distance(*(a+1)),trailw));
				tail.setOrigin(0,trailw/2);
				tail.setFillColor(sf::Color(rand() % 128,rand() % 64,rand() % 256,rand() % 128+128));
				tail.setPosition((*a).x(),(*a).y());
				tail.setRotation((-ang/Creature::PI)*180);
				trailw--;
				window->draw(tail);
			}
		}
		cToken->setFillColor(color);
		cToken->setRotation(creature->getAngleInDegrees()+90.f);
		cToken->setPosition(creature->getx(),creature->gety());
		ss.str("");
		ss << (int)creature->getFitness();
		text->setString(ss.str());
		sf::Vector2f tpos = cToken->getPosition();
		tpos.y -= cToken->getLocalBounds().height*2;
		text->setPosition(tpos);
		window->draw(*cToken);
		window->draw(*text);
		color.r -= 255/(creatures.size()*2);
		color.g -= 255/(creatures.size()*2);
		color.b -= 255/(creatures.size()*2);
	}
	for(auto& f:zoneTokens){
		window->draw(f);
	}
	window->display();
}
double frame_start, frame_end;
bool mainLoop() {
	frameclock->restart();
	bool runProgram = update();
	render();
	sf::Time frame_end = frameclock->getElapsedTime();
	sf::Time fps = sf::milliseconds(1000.f/60.f);
	if(frame_end<fps) {
		sf::sleep(fps-frame_end);
	}
	return runProgram;
}

int stringPopulation() {
	Population pop;
	pop.setTarget(-25);
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

