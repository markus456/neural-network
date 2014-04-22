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
#include <SFML\Graphics.hpp>
class Genome {
protected:
	std::string dna;
	std::vector<double> dna_values;
	double genomeFitness;
public:
	static const int DNA_LENGTH = 64;//Valuechain length
	static const int VALUE_BEGIN = 48;//First possible value
	static const int VALUE_END = 61;//Last possible value, inclusive
	static const int MUTATION_PROBABILITY = 25;//The probability of a value mutating to some random value
	static const int CROSSOVER = 70;
	Genome(std::vector<double>& init):
		dna("Genome")
	{
		genomeFitness = 1.f;
		dna_values = init;
		while (dna.length() < DNA_LENGTH) {
			dna += VALUE_BEGIN + rand() % (VALUE_END - VALUE_BEGIN + 1);
		}
	}
	Genome(unsigned int init):
		dna("Genome")
	{
		genomeFitness = 1.f;
		while (dna.length() < DNA_LENGTH) {
			dna += VALUE_BEGIN + rand() % (VALUE_END - VALUE_BEGIN + 1);
		}
	}
	Genome(std::string s = "") :
		dna(s)
	{
		genomeFitness = 1;
		while (dna.length() < DNA_LENGTH) {
			dna += VALUE_BEGIN + rand() % (VALUE_END - VALUE_BEGIN + 1);
		}
	}
	std::vector<double> getData() const
	{
		return dna_values;
	}
	std::string get()const
	{
		return dna;
	}
	virtual Genome operator+(const Genome& g)const
	{
		std::string child_dna;
		std::string pair = g.get();
		if(rand() % 2) {
			child_dna.append(dna.begin(), dna.begin() + (DNA_LENGTH / 2));
			child_dna.append(pair.begin() + (DNA_LENGTH / 2), pair.end());
		} else {
			child_dna.append(pair.begin(), pair.begin() + (DNA_LENGTH / 2));
			child_dna.append(dna.begin() + (DNA_LENGTH / 2), dna.end());
		}
		for (int i = 0; i < child_dna.length(); i++) {
			if (rand() % 100 < MUTATION_PROBABILITY) {
				char c = child_dna[i]+ (rand() % 3 - 1);
				c = c>VALUE_END? VALUE_BEGIN:c;
				c = c<VALUE_BEGIN?VALUE_END:c;
				child_dna.replace(rand() % child_dna.size(), 1, 1, c);
			}
		}
		return Genome(child_dna);
	}
	Genome& operator=(const Genome& g)
	{
		dna = g.dna;
		genomeFitness = g.genomeFitness;
		dna_values = g.getData();
		return *this;
	}
	double fitness()const
	{
		return genomeFitness;
	}
	void fitness(double f)
	{
		genomeFitness = f;
	}
	virtual ~Genome() {}
};
class Neuron {
private:
	int connections;
	std::vector<double> weights,in_values,prev_delta;
	double bias,threshold,out_value,signal_error,error_delta,momentum;

public:
	static const double ACTIVATION_RESPONSE;
	Neuron(int c = 0):
		connections(c), bias(-1.0), threshold(1.5), out_value(0.0), signal_error(0.0),error_delta(0.0),momentum(0.3)
	{
		std::mt19937 rnum_gen(std::chrono::system_clock::now().time_since_epoch().count());
		for(int i = 0; i<c+1; i++) {
			double rnd = (double)rnum_gen()/(double)rnum_gen.max() - (double)rnum_gen()/(double)rnum_gen.max();
			weights.push_back(rnd);
			in_values.push_back(0.0);
			prev_delta.push_back(0.0);
		}
	}
	void updateWeight(unsigned int k,double delta)
	{
		weights[k] += delta + momentum * prev_delta[k];
		prev_delta[k] = delta;
	}
	void setDelta(double d)
	{
		error_delta = d;
	}
	double getDelta()
	{
		return error_delta;
	}
	void setError(double err)
	{
		signal_error = err;
	}
	std::vector<double>& getInputValues()
	{
		return in_values;
	}
	double getError()
	{
		return signal_error;
	}
	unsigned int dnaSize()
	{
		return weights.size();
	}
	std::vector<double>& getWeights()
	{
		return weights;
	}
	void setWeights(std::vector<double> newvalues)
	{
		if(newvalues.size()==weights.size()) {
			weights = newvalues;
		}
	}
	void format(std::list<double>& new_w)
	{
		if(new_w.size()==weights.size())
		{
			for(int i = 0; i<weights.size(); i++) {
				weights[i] = new_w.front();
				new_w.pop_front();
			}
		}
	}
	double output(std::vector<double>& inputs)
	{
		double activation = 0;
		if(inputs.size()==connections) {
			for(int i = 0; i<weights.size()-1; i++) {
				activation += weights[i]*inputs[i];
				in_values[i] = inputs[i];
			}
			activation += weights[weights.size()-1]*bias;
		}
		//double th = tanh(activation) + 0.001*activation;
		double sigmoid = (1.f/(1+exp(-activation/ACTIVATION_RESPONSE)));
		out_value = sigmoid;
		return sigmoid;
	}
	double getOutputValue()
	{
		return out_value;
	}
};
class Membrane {
protected:
	std::vector<Neuron> neurons;
	int inputs,outputs;
public:
	Membrane(int nrns = 0, int cnns = 0)
	{
		for(int i = 0; i<nrns; i++) {
			neurons.push_back(Neuron(cnns));
		}
	}

	std::vector<Neuron>& getNeurons()
	{
		return neurons;
	}

	void format(std::list<double>& new_w)
	{
		for(int i = 0; i<neurons.size(); i++) {
			neurons[i].format(new_w);
		}
	}
	std::vector<double> process(std::vector<double>& inputs)
	{
		std::vector<double> output;
		for(auto& a:neurons) {
			output.push_back(a.output(inputs));
		}
		return output;
	}
	std::vector<double> getGenomes()
	{
		std::vector<double> genes;
		std::vector<double> input;
		for(int i = 0; i<neurons.size(); i++) {
			input = neurons[i].getWeights();
			genes.insert(genes.end(),input.begin(),input.end());
		}
		return genes;
	}
	unsigned int size()const
	{
		return neurons.size();
	}
	unsigned int layerDnaSize()
	{
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
	Brain(int in = 2,int out = 2,int lrs = 1,int nrn = 6):
		n_inputs(in),n_outputs(out)
	{
		layers.reserve(lrs);
		int neuron_count = nrn;
		int previousLayerSize = in;
		for(int i = 0; i<lrs; i++) {
			layers.push_back(Membrane(neuron_count,previousLayerSize));
			previousLayerSize = layers.back().size();
		}
		layers.push_back(Membrane(out,previousLayerSize));///output layer
	}
	unsigned int brainDnaSize()
	{
		unsigned int sz = 0;
		for(auto& a:layers) {
			sz += a.layerDnaSize();
		}
		return sz;
	}
	virtual ~Brain() {}
	void format(Genome new_dna)
	{
		auto data = new_dna.getData();
		std::list<double> nval;
		nval.assign(data.begin(),data.end());
		for(int i = 0; i<layers.size(); i++) {
			layers[i].format(nval);
		}
	}
	void format(std::list<double>& new_w)
	{
		for(int i = 0; i<layers.size(); i++) {
			layers[i].format(new_w);
		}
	}
	std::vector<double> update(std::vector<double>& input)
	{
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
	int nInputs()const
	{
		return n_inputs;
	}
	int nOutputs()const
	{
		return n_outputs;
	}

	std::vector<Membrane>& getMembranes()
	{
		return layers;
	}

	std::list<double> buildGenebase()
	{
		std::list<double> genebase;
		std::vector<double> tmp;
		for(int i = 0; i<layers.size(); i++) {
			tmp = layers[i].getGenomes();
			genebase.insert(genebase.end(),tmp.begin(),tmp.end());
		}
		return genebase;
	}
};

class Beholder:public Brain {
protected:
	unsigned int columns,rows,outputs;
	double learning_rate,base_rate,global_error;

public:
	///Number of pixels in inputs and number of output bits, number neuron layers and the neuron count for each layer.
	Beholder(int ccount, int rcount, unsigned int obits, int lrs, int nrns):Brain(ccount*rcount,obits,lrs,nrns), columns(ccount),rows(rcount),outputs(obits),learning_rate(0.1),base_rate(0.01) {}
	virtual void train(std::vector<double> inputs,char exp)
	{
		auto values = inputs;
		int xp = strtol(&exp,0,0);
		if(exp=='9'+1){
			xp = 10;
		}else if(exp=='9'+2){
			xp = 11;
		}
		int answer = analyze(values);
		std::vector<double> deltas;
		for(int i = 0; i<outputs; i++) {
			int var1 = xp&1?1:0;
			if(var1) {
				deltas.push_back(0.75);
			} else {
				deltas.push_back(0.25);
			}
			xp >>= 1;
		}
		///Output layer error
		int o_layer =getMembranes().size()-1;
		int z = 0;
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			double orig_out = getMembranes()[o_layer].getNeurons()[i].getOutputValue();
			double delta = deltas[z++] - orig_out;
			getMembranes()[o_layer].getNeurons()[i].setDelta(delta);
			double error = delta*orig_out*(1-orig_out);
			getMembranes()[o_layer].getNeurons()[i].setError(error);
		}
		///Hidder layer error
		for(int i = getMembranes().size()-2; i>-1; i--) {
			int x = 0;
			for(int j = getMembranes()[i].getNeurons().size()-1; j>-1; j--) {
				double sum = 0;
				for(int k = getMembranes()[i+1].getNeurons().size()-1; k>-1; k--) {
					sum += getMembranes()[i+1].getNeurons()[k].getError()*getMembranes()[i+1].getNeurons()[k].getWeights()[j];
				}
				double error = getMembranes()[i].getNeurons()[j].getOutputValue()*sum*(1-getMembranes()[i].getNeurons()[j].getOutputValue());
				getMembranes()[i].getNeurons()[j].setError(error);

			}
		}

		///Output layer correction
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			for(int j = 0; j<getMembranes()[o_layer].getNeurons()[i].getWeights().size(); j++) {
				double delta_w = learning_rate*getMembranes()[o_layer].getNeurons()[i].getInputValues()[j]*getMembranes()[o_layer].getNeurons()[i].getError();
				getMembranes()[o_layer].getNeurons()[i].updateWeight(j,delta_w);
			}
		}

		///Hidden layer correction
		for(int i = getMembranes().size()-2; i>-1; i--) {
			for(int j = getMembranes()[i].getNeurons().size()-1; j>-1; j--) {
				for(int k = 0; k<getMembranes()[i].getNeurons()[j].getWeights().size(); k++) {
					double delta_w = learning_rate*getMembranes()[i].getNeurons()[j].getInputValues()[k]*getMembranes()[i].getNeurons()[j].getError();
					getMembranes()[i].getNeurons()[j].updateWeight(k,delta_w);
				}
			}
		}
		global_error = 0;
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			global_error += pow(getMembranes()[o_layer].getNeurons()[i].getDelta(),2);
		}
		if(learning_rate>base_rate){
			learning_rate *= 0.9;
		}
	}
	virtual void train(std::vector<sf::RectangleShape>& inputs,char exp)
	{
		const char* expected = &exp;
		int xp = strtol(expected,0,0);
		int answer = analyze(inputs);
		std::vector<double> deltas;
		int combined = answer^xp;
		for(int i = 0; i<outputs; i++) {
			int var1 = xp&1?1:0;
			if(var1) {
				deltas.push_back(0.75);
			} else {
				deltas.push_back(0.25);
			}
			xp >>= 1;
		}
		///Output layer error
		int o_layer =getMembranes().size()-1;
		int z = 0;
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			double orig_out = getMembranes()[o_layer].getNeurons()[i].getOutputValue();
			double delta = deltas[z++] - orig_out;
			getMembranes()[o_layer].getNeurons()[i].setDelta(delta);
			double error = delta*orig_out*(1-orig_out);
			getMembranes()[o_layer].getNeurons()[i].setError(error);
		}
		///Hidder layer error
		for(int i = getMembranes().size()-2; i>-1; i--) {
			int x = 0;
			for(int j = getMembranes()[i].getNeurons().size()-1; j>-1; j--) {
				double sum = 0;
				for(int k = getMembranes()[i+1].getNeurons().size()-1; k>-1; k--) {
					sum += getMembranes()[i+1].getNeurons()[k].getError()*getMembranes()[i+1].getNeurons()[k].getWeights()[j];
				}
				double error = getMembranes()[i].getNeurons()[j].getOutputValue()*sum*(1-getMembranes()[i].getNeurons()[j].getOutputValue());
				getMembranes()[i].getNeurons()[j].setError(error);

			}
		}

		///Output layer correction
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			for(int j = 0; j<getMembranes()[o_layer].getNeurons()[i].getWeights().size(); j++) {
				double delta_w = (2*learning_rate)*getMembranes()[o_layer].getNeurons()[i].getInputValues()[j]*getMembranes()[o_layer].getNeurons()[i].getError();
				getMembranes()[o_layer].getNeurons()[i].updateWeight(j,delta_w);
			}
		}

		///Hidden layer correction
		for(int i = getMembranes().size()-2; i>-1; i--) {
			for(int j = getMembranes()[i].getNeurons().size()-1; j>-1; j--) {
				for(int k = 0; k<getMembranes()[i].getNeurons()[j].getWeights().size(); k++) {
					double delta_w = learning_rate*getMembranes()[i].getNeurons()[j].getInputValues()[k]*getMembranes()[i].getNeurons()[j].getError();
					getMembranes()[i].getNeurons()[j].updateWeight(k,delta_w);
				}
			}
		}
		global_error = 0;
		for(int i = getMembranes()[o_layer].getNeurons().size()-1; i>-1; i--) {
			global_error += pow(getMembranes()[o_layer].getNeurons()[i].getDelta(),2);
		}
		if(learning_rate>base_rate){
			learning_rate *= 0.9;
		}
	}
	void increaseLearningRate(){
		learning_rate *= 1.15;
	}
	void decreaseLearningRate(){
		if(learning_rate>base_rate){
			learning_rate = base_rate;
		}else{
			learning_rate *= 0.85;
		}
	}
	double getGlobalError()
	{
		return global_error;
	}
	virtual int analyze(std::vector<double>& inputs){
		if(inputs.size()==columns*rows) {
			auto outs = Brain::update(inputs);
			int result = 0;
			for(int i = 0; i<outs.size(); i++) {
				result <<= 1;
				result |= outs[i]>0.5?1:0;
			}
			return result;
		}
		return -1;
	}

	virtual int analyze(std::vector<sf::RectangleShape>& inputs)
	{
		if(inputs.size()==columns*rows) {
			std::vector<double> values;
			for(auto& a:inputs) {
				auto clr = a.getFillColor();
				if(clr==sf::Color::Black) {
					values.push_back(1.f);
				} else {
					values.push_back(0.f);
				}
			}

			auto outs = Brain::update(values);
			int result = 0;
			for(int i = 0; i<outs.size(); i++) {
				result <<= 1;
				result |= outs[i]>0.5?1:0;
			}
			return result;
		}
		return -1;
	}
};
template <typename T> class DataSet{
protected:
	std::vector<T> values;
	std::vector<unsigned int> indices;
public:
	void clear(){
		values.clear();
		indices.clear();
	}
	void reserve(unsigned int x, unsigned int y){
		values.reserve(x*y);
		indices.reserve(y);
	}
	unsigned int rows(){
		return indices.size();
	}
	unsigned int columns(unsigned int row){
		if(row<indices.size()-1){
			return indices[row+1] - indices[row];
		}else if(row<indices.size()){
			return indices[row] - indices[row-1];
		}
	}

	std::vector<T> at(unsigned int y){
		std::vector<T> row;
		if(y<indices.size()-1){
			row.assign(values.begin()+indices[y],values.begin()+indices[y+1]);
		}else if(y<indices.size()){
			row.assign(values.begin()+indices[y],values.end());
		}
		return row;
	}
	T& at(unsigned int x, unsigned int y){
		if(y<indices.size()-1){
			if(x<indices[y+1]-indices[y]){
				return values[indices[y]+x];
			}
		}else if(y<indices.size()){
			if(x+indices[y]<values.size()){
				return values[indices[y]+x];
			}
		}
		return T();
	}
	void push_row(std::vector<T> row){
		indices.push_back(values.size());
		values.insert(values.end(),row.begin(),row.end());
	}
	void push_back(T cell){
		values.push_back(cell);
	}
	void new_row(){
		indices.push_back(values.size());
	}
	void pop_back(){
		values.pop_back();
	}
	void pop_row(){
		values.erase(values.begin()+indices.back(),values.end());
		indices.pop_back();
	}
};

#endif
