#include "genetics.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <thread>
#include <sstream>
#include <SFML\Graphics.hpp>
#define GENERATIONS 1000
#define SCREEN_W 1200
#define SCREEN_H 800
#define MAX_TRAINING_ITERATIONS 5
#define X_RES 16
#define Y_RES 16
#define TILE_X_RES 32
#define TILE_Y_RES 32
#define IMAGE_DATA 1
#define IMAGE_DATA_W 4
#define IMAGE_DATA_H 3
bool update();
void render();
void train();
void log(std::string);
void exportBrains();
void importBrains();
std::vector<sf::RectangleShape> input_grid;
std::map<char,std::vector<std::string>> training_data;
std::map<char,DataSet<double>> training_set;
std::vector<char> found_sets;
sf::RenderWindow window;
sf::Clock beholder_speed;
sf::Image data;
sf::View input_camera;
sf::View log_camera;
sf::Font font;
sf::Vector2f logstart;
std::deque<sf::Text> log_output;
Beholder* beholder;
int fcnt,evolution_rate,increase_mutations, decrease_mutations,training_iters;
double speed, pavg, avg, current_mutation_rate,desired_error;
int counter, gen;
bool training;
unsigned int smallest_sample;
int main(int argc, char** argv)
{
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window.create(sf::VideoMode(SCREEN_W,SCREEN_H),"The Beholder",sf::Style::Default,settings);
	window.setFramerateLimit(60);


	input_camera.setViewport(sf::FloatRect(0.f,0.f,0.5f,1.f));
	log_camera.setViewport(sf::FloatRect(0.5f,0,0.5,1.f));


	training_iters = MAX_TRAINING_ITERATIONS;
	training = true;
	desired_error = 0.1;

	beholder = new Beholder(X_RES,Y_RES,4,2,(X_RES*Y_RES)/2);

	font.loadFromFile("revalia.ttf");

	if(IMAGE_DATA) {
		if(!data.loadFromFile("numbers.png")) {
			std::cout << "error loading image data!\n";
		}
		char target = '0';
		for(int ypos = 0; ypos<IMAGE_DATA_H; ypos++) {
			for(int xpos = 0; xpos<IMAGE_DATA_W; xpos++) {
				std::vector<double> symbol;
				for(int y = ypos*Y_RES; y<ypos*Y_RES+Y_RES; y++) {
					for(int x = xpos*X_RES; x<xpos*X_RES+X_RES; x++) {
						if(data.getPixel(x,y)==sf::Color::Black) {
							symbol.push_back(1.0);
						} else {
							symbol.push_back(0.0);
						}
					}
				}
				training_set[target].push_row(symbol);
				found_sets.push_back(target);
				target++;
			}
		}
	} else {
		std::stringstream ss;
		std::string fp;
		std::fstream infile;
		for(char start = '0'; start<'9'+1; start++) {
			ss.str("");
			ss << start << ".txt";
			fp = ss.str();
			infile.open(fp.c_str(),std::ios_base::in);
			smallest_sample = ~0;
			if(infile) {
				std::stringstream sstrm;
				sstrm << "training set for '" << start << "' found, reading to memory.\n";
				char buffer[X_RES*Y_RES+1];
				while(infile.getline(buffer,X_RES*Y_RES+1)) {
					training_data[start].push_back(buffer);
				}

				smallest_sample = training_data[start].size()<smallest_sample?training_data[start].size():smallest_sample;
				std::cout << "found " << training_data[start].size() << " sets of data.\n";
			} else {
				std::cout << "training set for '" << start << "' not found.\n";
			}
			infile.close();
		}
		std::cout << "smallest sample: " << smallest_sample << "\n";

		for(char start = '0'; start<'9'+1; start++) {
			if(training_data.count(start)>0) {
				found_sets.push_back(start);
			}
		}
		for(char ch:found_sets) {
			for(auto& line: training_data[ch]) {
				std::vector<double> row;
				for(auto line_iter = line.begin(); line_iter!=line.end(); line_iter++) {
					if(*line_iter=='1') {
						row.push_back(1.0);
					} else {
						row.push_back(0.0);
					}
				}
				std::random_shuffle(row.begin(),row.end());
				training_set[ch].push_row(row);
			}
		}
	}
	for(int y = 0; y<Y_RES; y++) {
		for(int x = 0; x<X_RES; x++) {
			input_grid.push_back(sf::RectangleShape(sf::Vector2f(32.f,32.f)));
			input_grid.back().setPosition(32.f*x,32.f*y);
			input_grid.back().setFillColor(sf::Color::White);
		}
	}
	sf::Vector2f sz;
	sz.x = input_grid.back().getPosition().x+input_grid.back().getLocalBounds().width-input_grid.front().getPosition().x;
	sz.y = input_grid.back().getPosition().y+input_grid.back().getLocalBounds().height-input_grid.front().getPosition().y;
	input_camera.setSize(sz);
	input_camera.setCenter(input_grid.front().getPosition().x+sz.x/2,input_grid.front().getPosition().y+sz.y/2);
	log_camera.setCenter(input_camera.getCenter());
	log_camera.move(0,sz.y);
	log_camera.setSize(sz);
	auto lpos = log_camera.getCenter();
	lpos.x -= log_camera.getSize().x/2;
	lpos.y -= log_camera.getSize().y/2;
	logstart = lpos;
	beholder_speed.restart();
	while(update())render();
	return 0;
}
void exportBrains(){
	if(beholder!=nullptr){
		std::fstream brains_out("saved.brains");
		std::stringstream ss;
		auto brains = beholder->buildGenebase();
		for(auto& a:brains){
			ss << a << '\n';
			auto str = ss.str();
			brains_out.write(str.c_str(),str.size());
		}
	}
}

void importBrains(){
	if(beholder!=nullptr){
		std::fstream brains_in("saved.brains");
		char buffer[256];
		std::vector<double> brains;
		while(brains_in.getline(buffer,256)){
			double w = strtod(buffer,0);
			brains.push_back(w);
		}
		Genome g(brains);
		beholder->format(g);
	}
}
void train()
{

	long iter = 0;
	int lc = 0;
	std::stringstream ss;
	std::vector<double> errors;
	errors.reserve(51);
	double error_avg = 9999.9;
	double prev_error = error_avg;
	log("training...\n");
	while(iter<training_iters) {
		for(char ch:found_sets) {
			beholder->train(training_set[ch].at(0),ch);
			errors.push_back(beholder->getGlobalError());
		}
		error_avg = 0.0;
		for(double d:errors) {
			error_avg += d;
		}
		error_avg *= 0.5;
		errors.clear();
		if(iter % 10 == 0) {
			ss.str("");
			ss << "Global error: " << error_avg << std::endl;
			log(ss.str());
			render();
		}
		iter++;
	}
	log("Training done!");

}
void log(std::string output){
	log_output.push_front(sf::Text(output,font,15));
	float ypos = logstart.y;
	for(auto& a = log_output.begin();a!=log_output.end();a++){
		(*a).setPosition(0,ypos);
		ypos += (*a).getLocalBounds().height;
		if(ypos - logstart.y>log_camera.getSize().y){
			log_output.erase(a,log_output.end());
			return;
		}
	}
}
bool update()
{
	sf::Event event;
	int ans = -1;
	std::stringstream ss;
	sf::String tmp;
	while(window.pollEvent(event)) {
		std::stringstream ss;
		switch (event.type) {
		case sf::Event::Closed:
			return false;
		case sf::Event::KeyPressed:
			switch(event.key.code){
			case sf::Keyboard::S:
				if(sf::Keyboard::isKeyPressed(sf::Keyboard::LControl)){
					exportBrains();
					log("Saved brains to file");
				}
				break;
			case sf::Keyboard::L:
				if(sf::Keyboard::isKeyPressed(sf::Keyboard::LControl)){
					importBrains();
					log("Loaded brains from file");
				}
				break;
			case sf::Keyboard::F1:
				for(auto& a:input_grid) {
					a.setFillColor(sf::Color::White);
				}
				break;
			case sf::Keyboard::F2:
				ans = beholder->analyze(input_grid);
				ss.str("");
				ss << "The Beholder sees a: " << ans;

				log(ss.str());
				break;
			case sf::Keyboard::F3:
				training = !training;

				ss.str("");
				ss << "Mode: ";
				if(training) {
					ss << "network training";
				} else {
					ss << "training set generation";
				}
				log(ss.str());
				break;
			case sf::Keyboard::F5:
				log("Training the beholder...");
				train();
				break;
			case sf::Keyboard::F9:

				ss.str("");
				training_iters += 5;
				ss << "Training iteration set to: " << training_iters << std::endl;
				log(ss.str());
				break;
			case sf::Keyboard::F10:
				if(training_iters>5){

					ss.str("");
					training_iters -= 5;
					ss << "Training iteration set to: " << training_iters << std::endl;
					log(ss.str());
				}

				break;
			}
			break;
		case sf::Event::TextEntered:
			if(IMAGE_DATA) {
				tmp = event.text.unicode;
				char c = tmp.toAnsiString()[0];
				if(c=='+') {
					c = '9'+1;
				} else if(c=='-') {
					c = '9'+2;
				}
				if(training_set.count(c)) {
					auto datarow = training_set[c].at(0);
					int i = 0;
					for(auto& a:input_grid) {
						if(datarow[i]==1.0) {
							a.setFillColor(sf::Color::Black);
						} else {
							a.setFillColor(sf::Color::White);
						}
						i++;
					}
				}
			} else {
				if(training) {
					tmp = event.text.unicode;
					beholder->train(input_grid,tmp.toAnsiString()[0]);
					ss << "The Beholder is trained with: " << tmp.toAnsiString() << " Global error: " << beholder->getGlobalError();
					log_output.push_front(sf::Text(sf::String(ss.str()),font,15));
				} else {
					tmp = event.text.unicode;
					ss << tmp.toAnsiString() << ".txt";
					std::string fpath = ss.str();
					std::fstream foutput(fpath.c_str(),std::ios_base::out|std::ios_base::app);
					ss.str("");
					for(auto& a:input_grid) {
						if(a.getFillColor()==sf::Color::Black) {
							ss << 1;
						} else {
							ss << 0;
						}
					}
					ss << '\n';
					auto str = ss.str();
					foutput.write(str.c_str(),str.size());
					ss.str("");
					ss << "Pattern saves as: " << tmp.toAnsiString() << " File is : " << fpath;
					log(ss.str());

					for(auto& a:input_grid) {
						a.setFillColor(sf::Color::White);
					}
				}
			}
			break;
		case sf::Event::MouseMoved:
			auto mp = window.mapPixelToCoords(sf::Mouse::getPosition(window),input_camera);
			if(sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
				for(auto& a:input_grid) {
					if(a.getGlobalBounds().contains(mp)) {
						a.setFillColor(sf::Color::Black);
						break;
					}
				}
			}
			if(sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
				for(auto& a:input_grid) {
					if(a.getGlobalBounds().contains(mp)) {
						a.setFillColor(sf::Color::White);
						break;
					}
				}
			}
			break;
		}
	}
	auto now = beholder_speed.getElapsedTime();
	if(now>sf::seconds(2)) {
		ans = beholder->analyze(input_grid);
		std::stringstream ss;
		ss << "The Beholder sees a: ";
		if(ans==10){
			ss << '+';
		}else if(ans==11){
			ss << '-';
		}else{
			ss << ans;
		}

		log(ss.str());
		beholder_speed.restart();
	}

	return true;
}
void render()
{
	window.clear();
	window.setView(input_camera);
	for(auto& a:input_grid) {
		window.draw(a);
	}
	window.setView(log_camera);
	for(auto& a:log_output) {
		window.draw(a);
	}
	window.display();
}
