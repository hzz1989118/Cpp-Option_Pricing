// This is a header file for EUoption Class
// Euro Option Class
// This class is based on Duffy's code
// parameters are encapulsed in this case, 
// whcih is easy for user to set and get its value

#ifndef EUoption_hpp
#define EUoption_hpp

#include "ExactOptionGlobal.hpp"
#include <boost/math/distributions/normal.hpp>
using namespace std;

class EUoption: public ExactOption
{
private:
	
	void copy(const EUoption& o2);
	double n(double x) const; //private functions has to be constant
	double N(double x) const;

	// 'Kernel' functions for option calculations
	double CallPrice(double U, const OptionPara& para) const;
	double PutPrice(double U, const OptionPara& para) const;
	double CallDelta(double U, const OptionPara& para) const;
	double PutDelta(double U, const OptionPara& para) const;
	double CallGamma(double U, const OptionPara& para) const;
	double PutGamma(double U, const OptionPara& para) const;


public:	// Public functions
	// No parametric constructor since member data is public
	EUoption();							// Default call option 
	EUoption(const EUoption& option2);	// Copy constructor
	EUoption (const OptionPara& para, string optionType);	// Specifier 
	virtual ~EUoption();	

	EUoption& operator = (const EUoption& option2);

	

	// Functions that calculate option price 
	double Price(double U) const; // U is current price for the underlying asset
	//vector and matrix pricer
	vector<double> Price(const vector<double>& U) const;
	virtual vector<vector<double> > Price(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const;
	

	// Sensitivities (Vector and Matrix Version)
	double Delta(double U) const;
	vector<double> Delta(const vector<double>& U) const;
	virtual vector<vector<double> > Delta(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const;
	
	double Gamma(double U) const;
	virtual vector<vector<double> > Gamma(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const;

	//Put-Call Parity
	double PriceUsingPutCallParity (double U) const; // Using PC Parity to calculate the call/put
	void CheckPutCallParity(double U) const; // Check the Parity

	//Exact Vs Numerical
	vector<double> ExactVSNumericalDelta(double U, const vector<double>& h) const; // For Delta Accuracy On each h
	vector<double> ExactVSNumericalGamma(double U, const vector<double>& h) const; // For Gamma Accuracy On each h
	//of course, this two can be extented to handle a vector of S (i.e., U)


};

#endif
