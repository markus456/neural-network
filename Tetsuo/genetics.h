#ifndef GENEGUARD
#define GENEGUARD
#include <string>
#include <cstdlib>
#include <utility>
#include <memory>
#include <vector>
#include <deque>
#include <map>
#include <list>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
class GenParams{
public:
	//Utility constants
	const static double PI;
	const static double PI_2;
	const static int TRAIL_LENGTH = 10;

	//Genome parameters
	const static int GENE_POOL_SIZE;
	static const int DNA_LENGTH = 64;//Valuechain length
	static const int VALUE_BEGIN = 48;//First possible value
	static const int VALUE_END = 61;//Last possible value, inclusive
	static const int MUTATION_PROBABILITY = 25;//The probability of a value mutating to some random value
	static const int CROSSOVER = 70;

	//Neural network parameters
	static const double ACTIVATION_RESPONSE;
};
class Point{
protected:
	double _x,_y,_w,_h;
public:
	Point(float x = 0, float y = 0, float w = 0, float h = 0):_x(x),_y(y),_w(w),_h(h){}
	Point& operator+(const Point& p){
		_x += p._x;
		_y += p._y;
		_w += p._w;
		_h += p._h;
		return *this;
	}
	Point& operator=(const Point& p){
		_x = p._x;
		_y = p._y;
		_w = p._w;
		_h = p._h;
		return *this;
	}
	Point& operator=(Point& p){
		_x = p._x;
		_y = p._y;
		_w = p._w;
		_h = p._h;
		return *this;
	}
	Point& operator=(Point&& p){
		_x = p._x;
		_y = p._y;
		_w = p._w;
		_h = p._h;
		return *this;
	}
	Point& operator=(const Point&& p){
		_x = p._x;
		_y = p._y;
		_w = p._w;
		_h = p._h;
		return *this;
	}
	Point& operator-(const Point& p){
		_x -= p._x;
		_y -= p._y;
		_w -= p._w;
		_h -= p._h;
		return *this;
	}
	double distance(const Point& p){
		return sqrt(pow(p._x-_x,2)+pow(p._y-_y,2));
	}
	//Radians
	double angRad(const Point& p){
		return atan2(-(p._y-_y),p._x-_x);
	}
	//Degrees
	double angDeg(const Point& p){
		return (atan2(-(p._y-_y),p._x-_x)/GenParams::PI)*180.0;
	}
	double x()const{return _x;}
	double y()const{return _y;}
	double w()const{return _w;}
	double h()const{return _h;}
};
class Genome {
protected:
	std::string dna;
	std::vector<double> dna_values;
	double genomeFitness;
public:
	
	Genome(std::vector<double>& init):dna("Genome") {
		genomeFitness = 1.f;
		dna_values = init;
		while (dna.length() < GenParams::DNA_LENGTH) {
			dna += GenParams::VALUE_BEGIN + rand() % (GenParams::VALUE_END - GenParams::VALUE_BEGIN + 1);
		}
	}
	Genome(const Genome& genome){
		genomeFitness = 1;
		dna_values = genome.dna_values;
		dna = genome.dna;
	}
	Genome(unsigned int init):dna("Genome") {
		genomeFitness = 1.f;
		while (dna.length() < GenParams::DNA_LENGTH) {
			dna += GenParams::VALUE_BEGIN + rand() % (GenParams::VALUE_END - GenParams::VALUE_BEGIN + 1);
		}
	}
	Genome(std::string s = "") :dna(s) {
		genomeFitness = 1;
		while (dna.length() < GenParams::DNA_LENGTH) {
			dna += GenParams::VALUE_BEGIN + rand() % (GenParams::VALUE_END - GenParams::VALUE_BEGIN + 1);
		}
	}
	std::vector<double> getData() const{
		return dna_values;
	}
	std::string get()const {
		return dna;
	}
	virtual Genome operator+(const Genome& g)const {
		std::string child_dna;
		std::string pair = g.get();
		if(rand() % 2) {
			child_dna.append(dna.begin(), dna.begin() + (GenParams::DNA_LENGTH / 2));
			child_dna.append(pair.begin() + (GenParams::DNA_LENGTH / 2), pair.end());
		} else {
			child_dna.append(pair.begin(), pair.begin() + (GenParams::DNA_LENGTH / 2));
			child_dna.append(dna.begin() + (GenParams::DNA_LENGTH / 2), dna.end());
		}
		for (int i = 0; i < child_dna.length(); i++) {
			if (rand() % 100 < GenParams::MUTATION_PROBABILITY) {
				char c = child_dna[i]+ (rand() % 3 - 1);
				c = c>GenParams::VALUE_END? GenParams::VALUE_BEGIN:c;
				c = c<GenParams::VALUE_BEGIN?GenParams::VALUE_END:c;
				child_dna.replace(rand() % child_dna.size(), 1, 1, c);
			}
		}
		return Genome(child_dna);
	}
	Genome& operator=(const Genome& g){
		dna = g.dna;
		genomeFitness = g.genomeFitness;
		dna_values = g.getData();
		return *this;
	}
	double fitness()const {
		return genomeFitness;
	}
	void fitness(double f) {
		genomeFitness = f;
	}
	virtual ~Genome() {}
};
class Population {
protected:
	std::vector<Genome> population;
	std::map<std::pair<double,double>, int> prop_val;
	double target_value;
	std::pair<std::string, double> closest;
	int POPULATION_SIZE;
public:

	Population(int p_size = 0) {
		POPULATION_SIZE = 0;
		if(p_size!=0) {
			POPULATION_SIZE = p_size;
			closest.second = 0;
			for (int i = 0; i < POPULATION_SIZE; i++) {
				Genome g;
				population.push_back(g);

			}
			closest.first = population.back().get();
			target_value = rand() % 100;
		}
	}
	void setTarget(double s) {
		target_value = s;
	}
	int getSize()const {
		return POPULATION_SIZE;
	}
	void clear() {
		population.clear();
		POPULATION_SIZE = 0;
	}
	void add(Genome g) {
		population.push_back(g);
		POPULATION_SIZE++;
	}
	virtual void evolve(std::vector<Genome>& nextgen, Genome& a, Genome& b) {
		std::string dna1 = a.get(),dna2 = b.get(),child1,child2;
		if(rand() % 100 < GenParams::CROSSOVER) {
			int rnd = rand() % GenParams::DNA_LENGTH;
			child1.append(dna1.begin(),dna1.begin()+rnd);
			child1.append(dna2.begin()+rnd,dna2.end());
			child2.append(dna2.begin(),dna2.begin()+rnd);
			child2.append(dna1.begin()+rnd,dna1.end());
		} else {
			child1 = dna1;
			child2 = dna2;
		}
		for (int i = 0; i < GenParams::DNA_LENGTH; i++) {
			if (rand() % 100 < GenParams::MUTATION_PROBABILITY) {
				char c = child1[i]+ (rand() % 3 - 1);
				c = c>GenParams::VALUE_END? GenParams::VALUE_BEGIN:c;
				c = c<GenParams::VALUE_BEGIN?GenParams::VALUE_END:c;
				child1.replace(i, 1, 1, c);
			}
			if (rand() % 100 < GenParams::MUTATION_PROBABILITY) {
				char c = child2[i]+ (rand() % 3 - 1);
				c = c>GenParams::VALUE_END? GenParams::VALUE_BEGIN:c;
				c = c<GenParams::VALUE_BEGIN?GenParams::VALUE_END:c;
				child2.replace(i, 1, 1, c);
			}
		}
		nextgen.push_back(Genome(child1));
		nextgen.push_back(Genome(child2));
	}
	virtual double evaluate(const Genome& g) {
		/*
		Characters to operators

		: = +
		; = -
		< = *
		= = /

		*/
		double diff = 0;
		int ptr = 0;
		std::list<int> values;
		std::list<char>operands;
		bool vnum = false,vopr = false;
		while(ptr<g.get().size()) {
			if(g.get()[ptr]-48<10) {
				if(vopr) {
					values.push_back(g.get()[ptr]-48);
					vopr = false;
				} else if(vnum) {
					int f = values.back();
					f = f*10 + (g.get()[ptr]-48);
					values.pop_back();
					values.push_back(f);
				} else {
					values.push_back(g.get()[ptr]-48);
				}
				vnum = true;
			} else if(vnum) {
				if(g.get()[ptr]-48>=10) {
					operands.push_back(g.get()[ptr]);
					vopr = true;
				}
			}
			ptr++;
		}
		if(values.size()>0) {
			diff += values.front();
			values.pop_front();
			while(values.size()>0&&operands.size()>0) {
				switch(operands.front()) {
				case ':':
					diff += values.front();
					values.pop_front();
					operands.pop_front();
					break;
				case ';':
					diff -= values.front();
					values.pop_front();
					operands.pop_front();
					break;
				case '<':
					diff = diff * values.front();
					values.pop_front();
					operands.pop_front();
					break;
				case '=':
					if(values.front()!=0) {
						diff = diff / values.front();
					}
					values.pop_front();
					operands.pop_front();
					break;
				}
			}
		}
		if(diff<0) {
			int d = -diff;
			return 1/(1+fabs(target_value+d));
		}
		return 1/(1+fabs(diff-target_value));
	}
	virtual void update() {
		double average = 0;
		std::vector<Genome> next_generation;
		for (auto& g : population) {
			g.fitness(evaluate(g));
			average += g.fitness();
		}
		double begin = 0;
		for(int i = 0; i<population.size(); i++) {
			double p = population[i].fitness()/average;
			p = p*100;
			prop_val.insert(std::make_pair(std::make_pair(begin,begin+p),i));
			begin += p;
		}
		while(next_generation.size()<POPULATION_SIZE) {
			int rnd, p1 = -1,p2 = -1;
			while(p1<0) {
				rnd = rand() % (int)(begin+1);
				for(auto& a:prop_val) {
					if(a.first.first<=rnd&&a.first.second>rnd) {
						p1 = a.second;
						break;
					}
				}
			}
			while(p2<0) {
				rnd = rand() % (int)(begin+1);
				for(auto& a:prop_val) {
					if(a.first.first<=rnd&&a.first.second>rnd) {
						if(a.second!=p1) {
							p2 = a.second;
						}
						break;
					}
				}
			}
			evolve(next_generation,population[p1],population[p2]);
		}

		prop_val.clear();
		population = next_generation;
	}

	virtual std::vector<Genome> getPopulation() {
		return population;
	}
	std::vector<Genome>& getPopulationRef() {
		return population;
	}
	double getTarget()const {
		return target_value;
	}
	virtual std::pair<std::string, double> fittestGenome() {
		double val = -1;
		for (auto& a : population) {
			if (evaluate(a) > val) {
				val = evaluate(a);
				closest.second = val;
				closest.first = a.get();
			}
		}
		return closest;
	}
};
class VectorPopulation:public Population {
protected:
	double random_factor;
public:
	VectorPopulation(int p_size = 0) {
		random_factor = 0.1;
		POPULATION_SIZE = 0;
		if(p_size!=0) {
			POPULATION_SIZE = p_size;
			closest.second = 0;
			for (int i = 0; i < POPULATION_SIZE; i++) {
				Genome g;
				population.push_back(g);

			}
			closest.first = population.back().get();

		}
		target_value = rand() % 100;
	}
	std::pair<std::string, double> fittestGenome() {
		return std::make_pair("null",0);
	}
	void setMutationScale(double d){
		random_factor = d;
	}
	void swapPopulation(std::vector<Genome>& genomes) {
		population.swap(genomes);
	}
	void evolve(std::vector<Genome>& nextgen, Genome& a, Genome& b) {
		std::vector<double> c1,c2,p1,p2;
		std::mt19937 rnum_gen(std::chrono::system_clock::now().time_since_epoch().count());
		p1 = a.getData();
		p2 = b.getData();
		if(p1.size()==0||p2.size()==0)return;
		if(rand() % 100 < GenParams::CROSSOVER) {
			for(int i = 0;i<p1.size();i++){
				int rnd = rand() % 100;
				float first = 100*(a.fitness()/(a.fitness()+b.fitness()));
				c1.push_back(rnd<first? p1.at(i):p2.at(i));
				c2.push_back(rnd<first? p2.at(i):p1.at(i));
			}
		} else {
			c1 = p1;
			c2 = p2;
		}
		if(rand() % 100< GenParams::MUTATION_PROBABILITY) {
			double rnd = (double)rnum_gen()/((double)rnum_gen.max()*0.99);
			rnd = rand() % 2? rnd:-rnd;
			c1[rand() %c1.size()] = rnd;
		}
		if(rand() % 100< GenParams::MUTATION_PROBABILITY) {
			double rnd = (double)rnum_gen()/((double)rnum_gen.max()*0.99);
			rnd = rand() % 2? rnd:-rnd;
			c2[rand() %c2.size()] = rnd;
		}
		for(int i = 0;i<c1.size();i++){
			if(rand() % 100< GenParams::MUTATION_PROBABILITY) {
				double rnd = (double)rnum_gen()/(double)rnum_gen.max();
				rnd = rand() % 2? rnd:-rnd;
				c1[i] += rnd;
			}
			if(rand() % 100< GenParams::MUTATION_PROBABILITY) {
				double rnd = (double)rnum_gen()/(double)rnum_gen.max();
				rnd = rand() % 2? rnd:-rnd;
				c2[i] += rnd;
			}
		}
		double af = (a.fitness()+b.fitness())/2;
		Genome g1(c1);
		g1.fitness(af);
		nextgen.push_back(g1);
	}
	double evaluate(const Genome& g) {
		return g.fitness();
	}
	void update() {
		double total = 0;
		std::vector<Genome> next_generation;
		std::map<std::string,double> propagate;
		for (auto& g : population) {
			propagate[g.get()] = evaluate(g);
			total += propagate[g.get()];
		}
		double valbegin = 0;
		for(int i = 0; i<population.size(); i++) {
			double p = population[i].fitness()/total;
			p = p*1000;
			prop_val.insert(std::make_pair(std::make_pair(valbegin,valbegin+p),i));
			valbegin += p;
		}
		if(POPULATION_SIZE==0&&population.size()>0) {
			POPULATION_SIZE = population.size();
		}
		std::sort(population.begin(),population.end(),[](const Genome& a,const Genome& b) {
			return a.fitness()<b.fitness();
		});
		while(next_generation.size()<POPULATION_SIZE/3){
			evolve(next_generation,population.front(),population.front());
		}
		while(next_generation.size()<POPULATION_SIZE) {
			int rnd, p1 = -1,p2 = -1;
			while(p1<0) {
				rnd = rand() % (int)(valbegin+1);
				for(auto& a:prop_val) {
					if(a.first.first<=rnd&&a.first.second>rnd) {
						p1 = a.second;
						break;
					}
				}
			}
			while(p2<0) {
				rnd = rand() % (int)(valbegin+1);
				for(auto& a:prop_val) {
					if(a.first.first<=rnd&&a.first.second>rnd) {
						if(a.second!=p1) {
							p2 = a.second;
						}
						break;
					}
				}
			}
			evolve(next_generation,population[p1],population[p2]);
		}

		prop_val.clear();
		population = next_generation;
	}
	std::vector<Genome> getPopulation() {
		return population;
	}
};
class Neuron {
private:
	int connections;
	std::vector<double> weights;
	double bias,threshold;

public:
	
	Neuron(int c = 0) {
		connections = c;
		threshold = 1.5;
		bias = -1.f;
		std::mt19937 rnum_gen(std::chrono::system_clock::now().time_since_epoch().count());
		for(int i = 0; i<c+1; i++) {
			double rnd = (double)rnum_gen()/((double)rnum_gen.max()*0.99);
			rnd = rand() % 2? rnd:-rnd;
			weights.push_back(rnd);
		}
	}
	unsigned int dnaSize() {
		return weights.size();
	}
	void format(std::list<double>& new_w) {
		for(int i = 0; i<weights.size(); i++) {
			weights[i] = new_w.front();
			new_w.pop_front();
		}
	}
	double output(std::vector<double>& inputs) {
		double activation = 0;
		if(inputs.size()==connections) {
			for(int i = 0; i<weights.size()-1; i++) {
				activation += weights[i]*inputs[i];
			}
			activation += weights[weights.size()-1]*bias;
		}
		double sigmoid = (1.f/(1+exp(-activation/GenParams::ACTIVATION_RESPONSE)));
		return sigmoid;
	}
	std::vector<double> getWeights() {
		return weights;
	}
};
class Membrane {
protected:
	std::vector<Neuron> neurons;
	int inputs,outputs;
public:
	Membrane(int nrns = 0, int cnns = 0) {
		for(int i = 0; i<nrns; i++) {
			neurons.push_back(Neuron(cnns));
		}
	}
	void format(std::list<double>& new_w) {
		for(int i = 0; i<neurons.size(); i++) {
			neurons[i].format(new_w);
		}
	}
	std::vector<double> process(std::vector<double>& inputs) {
		std::vector<double> output;
		for(auto& a:neurons) {
			output.push_back(a.output(inputs));
		}
		return output;
	}
	std::vector<double> getGenomes() {
		std::vector<double> genes;
		std::vector<double> input;
		for(int i = 0; i<neurons.size(); i++) {
			input = neurons[i].getWeights();
			genes.insert(genes.end(),input.begin(),input.end());
		}
		return genes;
	}
	unsigned int size()const {
		return neurons.size();
	}
	unsigned int layerDnaSize() {
		unsigned int sz = 0;
		for(auto& a:neurons) {
			sz += a.dnaSize();
		}
		return sz;
	}
};
class Brain {
protected:
	int n_inputs, n_outputs;
	std::vector<Membrane> layers;
public:
	Brain(int in = 2,int out = 2,int lrs = 1,int nrn = 6):n_inputs(in),n_outputs(out) {
		layers.reserve(lrs);
		int neuron_count = nrn;
		int previousLayerSize = in;
		for(int i = 0; i<lrs; i++) {
			layers.push_back(Membrane(neuron_count,previousLayerSize));
			previousLayerSize = layers.back().size();
		}
		layers.push_back(Membrane(out,previousLayerSize));///output layer
	}
	unsigned int brainDnaSize() {
		unsigned int sz = 0;
		for(auto& a:layers) {
			sz += a.layerDnaSize();
		}
		return sz;
	}
	virtual ~Brain() {}
	void format(Genome new_dna) {
		auto data = new_dna.getData();
		std::list<double> nval;
		nval.assign(data.begin(),data.end());
		for(int i = 0; i<layers.size(); i++) {
			layers[i].format(nval);
		}
	}
	void format(std::list<double>& new_w) {
		for(int i = 0; i<layers.size(); i++) {
			layers[i].format(new_w);
		}
	}
	std::vector<double> update(std::vector<double>& input) {
		std::vector<double> output;
		if(input.size()!=n_inputs)return output;
		for(int i = 0; i<layers.size(); i++) {
			if(i>0) {
				input.swap(output);
			}
			output = layers[i].process(input);
		}
		return output;
	}
	int nInputs()const {
		return n_inputs;
	}
	int nOutputs()const {
		return n_outputs;
	}
	std::list<double> buildGenebase() {
		std::list<double> genebase;
		std::vector<double> tmp;
		for(int i = 0; i<layers.size(); i++) {
			tmp = layers[i].getGenomes();
			genebase.insert(genebase.end(),tmp.begin(),tmp.end());
		}
		return genebase;
	}
};
class Creature {
protected:
	double x,y,vx,vy,target,angle,sizex,sizey,speed,tx,ty,l_leg,r_leg;
	Genome dna;
	VectorPopulation gene_evolution;
	Brain* b;
	int hitcount, n_inputs, n_outputs,neuron_count,layer_count;
	std::vector<double> output;
	std::deque<Point> trail;

public:
	Creature(int n_in = 4, int n_out = 2,int lrs = 1,int nrn_layer = 6) {
		hitcount = 0;
		sizex = sizey = 0;
		target = rand() % 180;
		target = (target/180)*GenParams::PI;
		speed = 2;
		l_leg = r_leg = 0.0;
		angle = target;
		vx = vy = 0;
		n_inputs = n_in;
		n_outputs = n_out;
		neuron_count = nrn_layer;
		layer_count = lrs;
		b = new Brain(n_in,n_out,lrs,nrn_layer);
		auto a = b->buildGenebase();
		std::vector<double> vals;
		vals.assign(a.begin(),a.end());
		dna = Genome(vals);
		setDNA(dna);
	}
	Creature(const Creature& c) {
		hitcount = 0;
		sizex = c.sizex;
		sizey = c.sizey;
		target = c.target;
		speed = c.speed;
		l_leg = r_leg = 0.0;
		angle = c.angle;
		vx = c.vx;
		vy = c.vy;
		x = c.x;
		y = c.y;
		n_inputs = c.n_inputs;
		n_outputs = c.n_outputs;
		neuron_count = c.neuron_count;
		layer_count = c.layer_count;
		b = new Brain(n_inputs,n_outputs,layer_count,neuron_count);
		dna = c.dna;
		setDNA(dna);
	}
	std::deque<Point>& getTrail(){return trail;}
	void foodAt(int x_,int y_) {
		target = atan2(-(y_-y),x_-x);
		tx = x_;
		ty = y_;
	}
	void areaSize(int x, int y) {
		sizex = x;
		sizey = y;
	}
	void setSpeed(double f) {
		speed = f;
	}
	void travelTo(double x_, double y_) {
		target = atan2(-(y_-y),x_-x);
	}
	double sigmoid(double f) {
		return 1.f/(1.f+exp(-f/GenParams::ACTIVATION_RESPONSE));
	}
	double normalize(double a, double b){
		return sqrt(pow(a,2)+pow(b,2));
	}
	void update() {

		if(trail.size()==0){
			trail.push_back(Point(x,y));
		}else if(trail.back().distance(Point(x,y))>10){
			trail.push_back(Point(x,y));
		}
		while(trail.size()>GenParams::TRAIL_LENGTH)trail.pop_front();
		std::vector<double> input;
		double ntx = tx-x, nty = -(ty-y), nx = vx, ny = vy;
		ntx = ntx/normalize(ntx,nty);
		nty = nty/normalize(ntx,nty);
		if(nx!=0||ny!=0){
			nx = nx/normalize(nx,ny);
			ny = ny/normalize(nx,ny);
		}
		/*if(target-GenParams::PI_2>angle){
			target = angle-GenParams::PI_2;
		}else if(target+GenParams::PI_2>angle){
			target = angle+GenParams::PI_2;
		}*/
		input.push_back(ntx);
		input.push_back(nty);
		input.push_back(nx);
		input.push_back(ny);
		output = b->update(input);
		l_leg = output[0];
		r_leg = output[1];
		double dir = l_leg-r_leg;
		dir = dir>GenParams::PI/8?GenParams::PI/8:dir;
		dir = dir<-GenParams::PI/8?-GenParams::PI/8:dir;
		double turnrate = (output[0]+output[1])/2.f;
		angle += dir;
		vx = cos(angle)*turnrate*speed;
		vy = sin(angle)*turnrate*speed;
		if(x+vx>0&&x+vx<sizex) {
			x += vx;
		}
		if(y-vy>0&&y-vy<sizey) {
			y -= vy;
		}
	}
	double targetX() {
		return tx;
	}
	double targetY() {
		return ty;
	}
	void increaseFitness(double f) {
		if(dna.fitness()+f>0) {
			dna.fitness(dna.fitness()+f);
		}
	}
	double getFitness() {
		return dna.fitness();
	}
	void setDNA(Genome g) {
		dna = g;
		b->format(dna);
	}
	Genome getDNA() {
		return dna;
	}
	std::vector<double> setInputs(std::vector<double>in) {
		if(in.size()==n_inputs) {
			output = b->update(in);
		} else {
			std::cout << "wrong number of inputs.\n";
		}
		return output;
	}
	void setPos(double x_, double y_) {
		x = x_;
		y = y_;
	}
	void setPos(Point p){
		x = p.x();
		y = p.y();
	}
	Point getPos(){
		return Point(x,y);
	}
	void setvx(double v) {
		vx = v;
	}
	void setvy(double v) {
		vy = v;
	}
	double getx()const {
		return x;
	}
	double gety()const {
		return y;
	}
	double getAngle() {
		return atan2(-vy,vx);
	}
	double getAngleInDegrees() {
		return (atan2(-vy,vx)/GenParams::PI)*180;
	}
	double getInvertedAngleInDegrees() {
		return (atan2(vy,vx)/GenParams::PI)*180;
	}
	double getvx()const {
		return vx;
	}
	double getvy()const {
		return vy;
	}
	std::vector<double> getOutputs()const {
		return output;
	}
	virtual ~Creature() {
		delete b;
	}
};
#endif
