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
sf::View message_camera;
sf::Font font;
sf::String kb_char;
sf::Text printer(sf::String(),font);
sf::Vector2f logstart;
std::deque<sf::Text> log_output;
sf::Text message;
Beholder* beholder;
Trainer* trainer;
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


	input_camera.setViewport(sf::FloatRect(0.f,0.f,0.5f,0.75f));
	message_camera.setViewport(sf::FloatRect(0.f,0.75f,0.5f,0.25f));
	log_camera.setViewport(sf::FloatRect(0.5f,0,0.5,1.f));


	training_iters = MAX_TRAINING_ITERATIONS;
	training = true;
	desired_error = 0.1;

	beholder = new Beholder(X_RES,Y_RES,7,2,(X_RES*Y_RES)/2);
	trainer = new Trainer("revalia.ttf");
	trainer->prepareData();
	trainer->assingStudent(beholder);
	font.loadFromFile("revalia.ttf");
	message = sf::Text("Draw Something!",font);
	message.setPosition(1000,1000);
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
	auto messagepos = message.getGlobalBounds();
	message_camera.setSize(sf::Vector2f(message.getLocalBounds().width,60));
	message_camera.setCenter(sf::Vector2f(messagepos.left+messagepos.width/2,messagepos.top+messagepos.height/2));
	beholder_speed.restart();
	while(update())render();
	return 0;
}
void exportBrains()
{
	if(beholder!=nullptr) {

		log("Exporting brains...");
		std::ofstream brains_out;
		brains_out.open("saved.txt");
		if(brains_out) {
			std::stringstream ss;
			std::vector<std::string> genestrings;
			auto brains = beholder->buildGenebase();
			std::cout <<  brains.size();
			for(auto& a:brains) {
				ss.str("");
				ss << a << '\n';
				brains_out << ss.str();
			}
			log( "Saved brains to file");
		} else {
			log("Error saving brains to file");
		}
		brains_out.close();
	}
}

void importBrains()
{
	if(beholder!=nullptr) {
		log("Importing brains...");
		std::ifstream brains_in("saved.txt");
		if(brains_in) {
			char buffer[256];
			std::list<double> brains;
			while(brains_in.getline(buffer,256)) {
				double w = strtod(buffer,0);
				brains.push_back(w);
			}
			beholder->format(brains);
			log("Loaded brains from file");
		} else {
			log("Error loading brains from file");
		}
		brains_in.close();
	}
}
void train()
{
	long iter = 0;
	char hardest = '0';
	double error_avg = 9999.9,prev_error = error_avg,hardest_error = 0.0;
	std::stringstream ss;
	std::vector<double> errors;
	errors.reserve(found_sets.size());
	message.setString("Training the Beholder...\n");
	while(iter<training_iters) {
		hardest_error = 0.0;
		for(char ch:found_sets) {
			beholder->train(training_set[ch].at(0),ch);
			errors.push_back(beholder->getGlobalError());
			if(errors.back()>hardest_error) {
				hardest_error = errors.back();
				hardest = ch;
			}
		}
		beholder->train(training_set[hardest].at(0),hardest);
		errors.push_back(beholder->getGlobalError());
		beholder->correct();
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
	message.setString("Training done!");

}
void log(std::string output)
{
	log_output.push_front(sf::Text(output,font,15));
	float ypos = logstart.y;
	for(auto a = log_output.begin(); a!=log_output.end(); a++) {
		(*a).setPosition(0,ypos);
		ypos += (*a).getLocalBounds().height;
		if(ypos - logstart.y>log_camera.getSize().y) {
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
		std::vector<double> board;
		switch (event.type) {
		case sf::Event::Closed:
			return false;
		case sf::Event::KeyPressed:
			switch(event.key.code) {
			case sf::Keyboard::S:
				if(sf::Keyboard::isKeyPressed(sf::Keyboard::LControl)) {
					exportBrains();

				}
				break;
			case sf::Keyboard::L:
				if(sf::Keyboard::isKeyPressed(sf::Keyboard::LControl)) {
					importBrains();

				}
				break;
			case sf::Keyboard::F1:
				for(auto& a:input_grid) {
					a.setFillColor(sf::Color::White);
				}
				break;
			case sf::Keyboard::F2:
				for(auto& a:input_grid) {
					if(a.getFillColor()== sf::Color::Black) {
						board.push_back(1.0);
					} else {
						board.push_back(0.0);
					}
				}
				message.setString(sf::String((char)beholder->analyze(board)));

				if(message.getLocalBounds().width>message_camera.getSize().x) {
					auto messagepos = message.getGlobalBounds();
					message_camera.setSize(sf::Vector2f(message.getLocalBounds().width,60));
					message_camera.setCenter(sf::Vector2f(messagepos.left+messagepos.width/2,messagepos.top+messagepos.height/2));
				}
				break;
			case sf::Keyboard::F3:
				training = !training;
				ss << "Mode: ";
				if(training) {
					ss << "network training";
				} else {
					ss << "training set generation";
				}
				log(ss.str());
				break;
			case sf::Keyboard::F5:
				for(int i = 0;i<25;i++)trainer->train();
				ss <<  "Error: " << trainer->globalError();
				log("Beholder trained 25 times");
				log(ss.str());
				break;
			case sf::Keyboard::F9:

				ss.str("");
				training_iters += 5;
				ss << "Training iteration set to: " << training_iters << std::endl;
				log(ss.str());
				break;
			case sf::Keyboard::F10:
				if(training_iters>5) {
					training_iters -= 5;
					ss << "Training iteration set to: " << training_iters << std::endl;
					log(ss.str());
				}

				break;
			}
			break;
		case sf::Event::TextEntered:
			kb_char = sf::String(event.text.unicode);


			break;
		}
	}

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

	window.setView(message_camera);
	window.draw(message);

	printer.setString(kb_char);
	printer.setColor(sf::Color::Black);
	printer.setOrigin(printer.getLocalBounds().width/2.f,printer.getLocalBounds().height/2.f);
	printer.setPosition(16.f,16.f);
	window.draw(printer);

	window.display();
}
