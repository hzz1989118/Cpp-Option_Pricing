// This is a header file for Perpetual American Option Class
// This class is based on Duffy's code
// parameters are encapulsed in this case, 
// whcih is easy for user to set and get its value

#ifndef PAMoption_hpp
#define PAMoption_hpp

#include "ExactOptionGlobal.hpp" // Get OptionPara Struct
// #include <boost/math/distributions/normal.hpp>
using namespace std;

class PAMoption: public ExactOption
{
private:

	void copy(const PAMoption& o2);

	// 'Kernel' functions for option calculations
	double CallPrice(double U, const OptionPara& para) const;
	double PutPrice(double U, const OptionPara& para) const;


public:	// Public functions
	// No parametric constructor since member data is public
	PAMoption();							// Default call option 
	PAMoption(const PAMoption& option2);	// Copy constructor
	PAMoption (const OptionPara& para, string optionType);	// Specifier by OptionPara and optionType
	virtual ~PAMoption();	
	// Assign Operator
	PAMoption& operator = (const PAMoption& option2);

	// Functions that calculate option price 
	double Price(double U) const; // U is current price for the underlying asset
	//vector and matrix pricer
	vector<double> Price(const vector<double>& U) const;
	virtual vector<vector<double> > Price(const vector<double>& U, const vector<vector<double> > & ParaMatrix) const;

};

#endif