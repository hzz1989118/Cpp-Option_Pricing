#ifndef ExactOptionGlobal_hpp
#define ExactOptionGlobal_hpp

#include<vector>
#include<string>
#include<random>
#include<iostream>
#include<boost/random.hpp>
#include <iomanip>

using namespace std;



// Data Encapsulate Struct
struct OptionPara
{
	double T;		// Expiry date
	double K;		// Strike price
	double sig;		// Volatility
	double r;		// Interest rate
	double b;		// Cost of carry
};



// BASE CLASS ExactOption

class ExactOption{
	private:

	OptionPara m_para;
	string optType;		

	void init(){	// Initialise all default values
		// Default values
		m_para.r = 0.05;
		m_para.sig = 0.2;
		m_para.K = 110.0;
		m_para.T = 0.5;
		m_para.b = m_para.r;			// Black and Scholes stock option model (1973)
		optType = "C";		// European Call Option (this is the default type)
	}

	public:	// Public functions

	ExactOption (){
		init();
	}
	ExactOption (const OptionPara& para, string optionType = "C"){
		m_para = para;
		optType = optionType;
		if (optType == "c")
			optType = "C";
	}	// Specifier 
	virtual ~ExactOption(){}	

	virtual vector<vector<double> > Price(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const{
		vector<vector<double> > out;
		return out;
	}	
	virtual vector<vector<double> > Delta(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const{
		vector<vector<double> > out;
		return out;
	}
	virtual vector<vector<double> > Gamma(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const{
		vector<vector<double> > out;
		return out;
	}	

	//Getter
	OptionPara GetPara() const{
		return(m_para);
	}

	string GetType() const{
		return(optType);
	}

	// Modifier functions
	void toggle(){
	// Change option type (C/P, P/C)

	if (optType == "C")
		optType = "P";
	else
		optType = "C";
	}
			// Change option type (C/P, P/C)
	void SetPara(const OptionPara& para){ 
	//set parametrs by struct OptionPara
		m_para = para;
	}
	void SetType(const string& optionType){
		optType = optionType;
	}
};

//Interfaces

//Printer
inline void Print1d(const vector<double>& v) // Print a double vector
{
	for(int i = 0; i < v.size(); i++){
		cout << v[i] << " ";
	}
}

inline void Print2d(const vector<vector<double> >& v) // Print a matrix
{	
	int d1 = v[0].size();
	int d2 = v.size();
	for(int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			cout << std::setprecision(3) << v[j][i] << " ";
		}
		cout << endl;
	}
	
}

// Mesh Generator
inline vector<double> getMesh(double start, double end, double h){
	// Global Function of getting Meshes by size
	vector<double> out;
	while (start <= end){
		out.push_back(start);
		start += h;
	}
	return out;
}

static std::random_device rd;
static boost::mt19937 gen(rd());

inline vector<double> getMesh(double start, double end, int n_rand){
	// Global Function of getting Meshes by unif random numbers
	vector<double> out;
	int i = 0;
	boost::random::uniform_real_distribution<> runif(start, end);
	while ( i < n_rand){
		out.push_back(runif(gen));
		i += 1;
	}
	return out;
}

inline vector<double> getMesh(){
	// Global Function for getting Meshes by console input
	std::vector<double> out;
	double input;
	cout << "Please give a number and type in anything else to stop: " << endl;
	while (cin >> input){
		out.push_back(input);
		cin.clear();
		cout << "Please give a number and type in anything else to stop: " << endl;
	}
	return out;
}

// Pricer Functions


inline vector<vector<double> > OptionMatrixPricer (const ExactOption& option, const vector<double>& U, const vector<vector<double> > & ParaMatrix){
	return option.Price(U, ParaMatrix);
}

inline vector<vector<double> > OptionMatrixDelta (const ExactOption& option, const vector<double>& U, const vector<vector<double> > & ParaMatrix){
	return option.Delta(U, ParaMatrix);
}

inline vector<vector<double> > OptionMatrixGamma (const ExactOption& option, const vector<double>& U, const vector<vector<double> > & ParaMatrix){
	return option.Gamma(U, ParaMatrix);
}




#endif