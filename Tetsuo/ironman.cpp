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
bool mainLoop();
int stringPopulation();
std::vector<std::shared_ptr<Creature>> creatures;
VectorPopulation* vpop;
std::vector<Point> food;
sf::Shape* cToken;
sf::Shape* fToken;
sf::RenderWindow* window;
sf::Clock* frameclock;
sf::Text* text;
int fcnt,evolution_rate,increase_mutations, decrease_mutations;
double speed, pavg, avg, current_mutation_rate,bestavg;
int counter, gen, decline_counter;
std::vector<Genome> best_so_far;
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
			std::sort(food.begin(),food.end(),[&](const Point& a,const Point& b) {
				return pow(a.x()-creature->getx(),2)+pow(a.y()-creature->gety(),2)<pow(b.x()-creature->getx(),2)+pow(b.y()-creature->gety(),2);
			});
			creature->foodAt(food.front().x(),food.front().y());
			creature->update();
		}
		for(auto& creature : creatures) {
			std::sort(food.begin(),food.end(),[&](const Point& a,const Point& b) {
				return pow(a.x()-creature->getx(),2)+pow(a.y()-creature->gety(),2)<pow(b.x()-creature->getx(),2)+pow(b.y()-creature->gety(),2);
			});
			if(sqrt(pow(food.front().x()-creature->getx(),2)+pow(food.front().y()-creature->gety(),2))<25) {
				creature->increaseFitness(1.f);
				food.erase(food.begin());
				//fcounter = 60*evolution_rate;
			}
		}
		while(food.size()<FOOD_CAP) {
			food.push_back(Point(20 + rand() % (SCREEN_W-40), 10 + rand() % (SCREEN_H-20)));
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
	gen = 1;
	evolution_rate = 10;
	counter = 300;
	bestavg = 0;
	decline_counter = 0;
	srand(time(0));
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window = new sf::RenderWindow(sf::VideoMode(SCREEN_W,SCREEN_H),"Evolution 0.3",sf::Style::Default,settings);
	sf::CircleShape* tmp = new sf::CircleShape(15.f,3);
	tmp->setOutlineColor(sf::Color::White);
	tmp->setFillColor(sf::Color::Red);
	tmp->setOutlineThickness(2.f);
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
		creatures.push_back(std::shared_ptr<Creature>(new Creature(4,2,1,6)));
		creatures.back()->setPos(10 + rand() % (SCREEN_W-20), 10 + rand() % (SCREEN_H-20));
		creatures.back()->areaSize(SCREEN_W,SCREEN_H);
		creatures.back()->setSpeed(4);
	}
	for(int i = 0; i<FOOD_CAP; i++) {
		food.push_back(Point(20 + rand() % (SCREEN_W-40), 10 + rand() % (SCREEN_H-20)));
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
				for(int i = 0;i<10;i++){
					updateCreatures();
				}
				break;
			case sf::Keyboard::F2:
				for(int i = 0;i<100;i++){
					updateCreatures();
				}
				break;
			case sf::Keyboard::F3:
				//fastForward(10);
				break;
			case sf::Keyboard::F4:
				//fastForward(50);
				break;
			}
			break;
		case sf::Event::MouseButtonPressed:

			Point p(event.mouseButton.x,event.mouseButton.y);
			auto a = creatures.begin();
			while(a!=creatures.end()){
				if(p.distance(Point((*a)->getx(),(*a)->gety()))<15){
					p = (*a)->getPos();
					a = creatures.erase(a);
					std::sort(creatures.begin(),creatures.end(),[](const std::shared_ptr<Creature>& a,const std::shared_ptr<Creature>& b){
						return a->getFitness()>b->getFitness();
					});
					if(event.mouseButton.button==sf::Mouse::Left){
						creatures.push_back(std::shared_ptr<Creature>(new Creature(*creatures.front())));
					}else if(event.mouseButton.button==sf::Mouse::Right){
						creatures.push_back(std::shared_ptr<Creature>(new Creature));
						creatures.back()->areaSize(SCREEN_W,SCREEN_H);
						creatures.back()->setSpeed(4);
					}
					creatures.back()->setPos(p);
					break;
				}else{
					a++;
				}
			}
			break;
		}
	}
	double ftest = 0;
	avg = 0;
	for(auto& creature: creatures) {
		std::sort(food.begin(),food.end(),[&](const Point& a,const Point& b) {
			return pow(a.x()-creature->getx(),2)+pow(a.y()-creature->gety(),2)<pow(b.x()-creature->getx(),2)+pow(b.y()-creature->gety(),2);
		});
		creature->foodAt(food.front().x(),food.front().y());
		creature->update();
		avg += creature->getFitness();
		if(creature->getFitness()>ftest){
			ftest = creature->getFitness();
		}
	}
	pavg = avg;
	avg = avg/creatures.size();
	if(avg>bestavg){
		best_so_far = vpop->getPopulationRef();
	}
	for(auto& creature : creatures) {
		std::sort(food.begin(),food.end(),[&](const Point& a,const Point& b) {
			return pow(a.x()-creature->getx(),2)+pow(a.y()-creature->gety(),2)<pow(b.x()-creature->getx(),2)+pow(b.y()-creature->gety(),2);
		});
		if(sqrt(pow(food.front().x()-creature->getx(),2)+pow(food.front().y()-creature->gety(),2))<25) {
			fcnt = 0;
			creature->increaseFitness(1.f);
			counter = 60*evolution_rate;
			food.erase(food.begin());
		}
	}

	while(food.size()<FOOD_CAP) {
		food.push_back(Point(20 + rand() % (SCREEN_W-40), 10 + rand() % (SCREEN_H-20)));
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
void render() {
	std::stringstream ss;
	std::deque<Point> trails, vertices;
	window->clear();


	std::sort(creatures.begin(),creatures.end(),[](const std::shared_ptr<Creature>& a,const std::shared_ptr<Creature>& b){
		return a->getFitness()>b->getFitness();
	});
	sf::Color color(sf::Color::Red);
	sf::ConvexShape tail;
	tail.setFillColor(sf::Color(128,0,0,128));
	for(auto& creature:creatures){
		trails = creature->getTrail();
		int pcount = 0;
		double trailw = 2;
		for(int a = 0;a<trails.size()-1;a++){
			double ang = trails[a].angRad(trails[a+1]);
			vertices.push_back(Point(trails[a].x()+cos(ang-GenParams::PI_2)*trailw/2,trails[a].y()-sin(ang-GenParams::PI_2)*trailw/2));
			vertices.push_back(Point(trails[a].x()+cos(ang+GenParams::PI_2)*trailw/2,trails[a].y()-sin(ang+GenParams::PI_2)*trailw/2));
			trailw++;
			vertices.push_back(Point(trails[a+1].x()+cos(ang+GenParams::PI_2)*trailw/2,trails[a+1].y()-sin(ang+GenParams::PI_2)*trailw/2));
			vertices.push_back(Point(trails[a+1].x()+cos(ang-GenParams::PI_2)*trailw/2,trails[a+1].y()-sin(ang-GenParams::PI_2)*trailw/2));
			int p = 0;
			tail.setPointCount(vertices.size());
			for(auto& a:vertices){
				tail.setPoint(p,sf::Vector2f(a.x(),a.y()));
				p++;
			}
			window->draw(tail);

			vertices.clear();
		}

		cToken->setFillColor(color);
		cToken->setRotation(creature->getAngleInDegrees()+90.f);
		cToken->setPosition(creature->getx(),creature->gety());
		ss.str("");
		ss << creature->getFitness();
		text->setString(ss.str());
		sf::Vector2f tpos = cToken->getPosition();
		tpos.y -= cToken->getLocalBounds().height*2;
		text->setPosition(tpos);
		window->draw(*cToken);
		window->draw(*text);
	}
	for(auto& f:food){
		fToken->setPosition(f.x(),f.y());
		window->draw(*fToken);
	}
	ss.str("");
	ss << "Generation: " << gen;
	text->setPosition(10,10);
	text->setString(ss.str());
	window->draw(*text);
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

