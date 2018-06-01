// EUoption.cpp
//
// Constructors and Member Functions of EUoption Class


#include "EUoption.hpp"
#include <cmath>
#include <iostream>

//Gaussian Pdf and Cdf

double EUoption::n(double x) const	
{ 
	boost::math::normal_distribution<double> EUoptionStandardNormal;
	return(pdf(EUoptionStandardNormal, x)); // I use boost package

};

double EUoption::N(double x) const
{ 	
	boost::math::normal_distribution<double> EUoptionStandardNormal;
	return(cdf(EUoptionStandardNormal, x)); // I use boost package
};



// Kernel Functions (Haug)
double EUoption::CallPrice(double U, const OptionPara& para) const
{
	double _T = para.T;
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double tmp = _sig * sqrt(_T);
	double d1 = ( log(U/_K) + (_b + (_sig*_sig)*0.5 ) * _T )/ tmp;
	double d2 = d1 - tmp;
	return (U * exp((_b-_r)*_T) * N(d1)) - (_K * exp(-_r * _T)* N(d2));

}

double EUoption::PutPrice(double U, const OptionPara& para) const
{

	double _T = para.T;
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double tmp = _sig * sqrt(_T);
	double d1 = ( log(U/_K) + (_b + (_sig*_sig)*0.5 ) * _T )/ tmp;
	double d2 = d1 - tmp;;
	return (_K * exp(-_r * _T)* N(-d2)) - (U * exp((_b-_r)*_T) * N(-d1));

}

double EUoption::CallDelta(double U, const OptionPara& para) const
{	
	double _T = para.T;
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double tmp = _sig * sqrt(_T);
	double d1 = ( log(U/_K) + (_b + (_sig*_sig)*0.5 ) * _T )/ tmp;
	return exp((_b-_r)*_T) * N(d1);
}

double EUoption::PutDelta(double U, const OptionPara& para) const
{
	double _T = para.T;
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double tmp = _sig * sqrt(_T);
	double d1 = ( log(U/_K) + (_b + (_sig*_sig)*0.5 ) * _T )/ tmp;

	return exp((_b-_r)*_T) * (N(d1) - 1.0);
}

double EUoption::CallGamma(double U, const OptionPara& para) const
{	
	double _T = para.T;
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double tmp = _sig * sqrt(_T);
	double d1 = ( log(U/_K) + (_b + (_sig*_sig)*0.5 ) * _T )/ tmp;
	double d2 = d1 - tmp;
	return ( n(d1) * exp((_b-_r)*_T) ) / (U * tmp);
}

double EUoption::PutGamma(double U, const OptionPara& para) const
{
	return CallGamma(U, para);
}


/////////////////////////////////////////////////////////////////////////////////////

void EUoption::copy( const EUoption& o2)
{

	this->SetPara(o2.GetPara());
	this->SetType(o2.GetType());
	
}

EUoption::EUoption(){
	ExactOption();
}

EUoption::EUoption(const OptionPara& para, string optionType){
	//specific constructor
	this->SetPara(para);
	this->SetType(optionType);
	if (this->GetType() == "c")
		this->SetType("C");
}

EUoption::EUoption(const EUoption& o2)
{ // Copy constructor
	copy(o2);
}


EUoption::~EUoption()
{

}


EUoption& EUoption::operator = (const EUoption& option2)
{

	if (this == &option2) return *this;

	copy (option2);

	return *this;
}

// Functions that calculate option price and sensitivities
double EUoption::Price(double U) const
{
	if (this->GetType() == "C")
	{	
		//cout << "calling call\n";
		return CallPrice(U, this->GetPara()); //uing member parameters
	}
	else
	{
		//cout << "calling put\n";
		return PutPrice(U, this->GetPara());
	}
}

//Overload Functions for Pricer
// Vector Pricer
vector<double> EUoption::Price(const vector<double>& U) const
{	
	int size_v = U.size();
	vector<double> out(size_v);
	for (int i = 0; i <  size_v; i++){
		out[i] = Price(U[i]);
	}
	return(out);
}

// Matrix Pricer
vector<vector<double> > EUoption::Price(const vector<double>& U, const vector<vector<double> >& ParaMatrix) const
{	
	int row_size = ParaMatrix.size(); // row size of output matrix
	int num_para = ParaMatrix[0].size(); // col num of input matrix
	int col_size = U.size(); // col size of output matrix

	if (num_para != 5)
	{	// To see if the input matrix is valid
		cout << "The parameter matrix requires 5 columns!" << endl;
		cout << "Colums are T, K, sig, r, b" << endl;
		exit(1);
	}

	vector<vector<double> > out(row_size, vector<double> (col_size));
	for (int i = 0; i < row_size; i++)
	{	
		// the column order does matter
		OptionPara TmpPara;
		TmpPara.T = ParaMatrix[i][0];
		TmpPara.K = ParaMatrix[i][1];
		TmpPara.sig = ParaMatrix[i][2];
		TmpPara.r = ParaMatrix[i][3];
		TmpPara.b = ParaMatrix[i][4];
		for(int j=0; j < col_size; j++)
		{
			if (this->GetType() == "C")
			{
				out[i][j] = CallPrice(U[j], TmpPara);
			} else {
				out[i][j] = PutPrice(U[j], TmpPara);
			}
		}
		
	}
	return(out);
}	

//Overload Functions for Delta
// Vector Delta

double EUoption::Delta(double U) const 
{
	if (this->GetType() == "C")
		return CallDelta(U, this->GetPara());
	else
		return PutDelta(U, this->GetPara());

}

vector<double> EUoption::Delta(const vector<double>& U) const
{	
	int size_v = U.size();
	vector<double> out(size_v);
	for (int i = 0; i <  size_v; i++){
		out[i] = Delta(U[i]);
	}
	return(out);
}

vector<vector<double> > EUoption::Delta(const vector<double>& U, const vector<vector<double> >& ParaMatrix) const
{	
	int row_size = ParaMatrix.size(); // row size of output matrix
	int num_para = ParaMatrix[0].size(); // col num of input matrix
	int col_size = U.size(); // col size of output matrix

	if (num_para != 5)
	{	// To see if the input matrix is valid
		cout << "The parameter matrix requires 5 columns!" << endl;
		cout << "Colums are T, K, sig, r, b" << endl;
		exit(1);
	}

	vector<vector<double> > out(row_size, vector<double> (col_size));
	for (int i = 0; i < row_size; i++)
	{	
		// the column order does matter
		OptionPara TmpPara;
		TmpPara.T = ParaMatrix[i][0];
		TmpPara.K = ParaMatrix[i][1];
		TmpPara.sig = ParaMatrix[i][2];
		TmpPara.r = ParaMatrix[i][3];
		TmpPara.b = ParaMatrix[i][4];
		for(int j=0; j < col_size; j++)
		{
			if (this->GetType() == "C")
			{
				out[i][j] = CallDelta(U[j], TmpPara);
			} else {
				out[i][j] = PutGamma(U[j], TmpPara);
			}
		}
		
	}
	return(out);
}

//Overload Functions for Gamma


double EUoption::Gamma(double U) const
{
	if (this->GetType() == "C")
		return CallGamma(U, this->GetPara());
	else
		return PutGamma(U, this->GetPara());
}

// Matrix Delta

vector<vector<double> > EUoption::Gamma(const vector<double>& U, const vector<vector<double> >& ParaMatrix) const
{	
	int row_size = ParaMatrix.size(); // row size of output matrix
	int num_para = ParaMatrix[0].size(); // col num of input matrix
	int col_size = U.size(); // col size of output matrix

	if (num_para != 5)
	{	// To see if the input matrix is valid
		cout << "The parameter matrix requires 5 columns!" << endl;
		cout << "Colums are T, K, sig, r, b" << endl;
		exit(1);
	}

	vector<vector<double> > out(row_size, vector<double> (col_size));
	for (int i = 0; i < row_size; i++)
	{	
		// the column order does matter
		OptionPara TmpPara;
		TmpPara.T = ParaMatrix[i][0];
		TmpPara.K = ParaMatrix[i][1];
		TmpPara.sig = ParaMatrix[i][2];
		TmpPara.r = ParaMatrix[i][3];
		TmpPara.b = ParaMatrix[i][4];
		for(int j=0; j < col_size; j++)
		{
			if (this->GetType() == "C")
			{
				out[i][j] = CallGamma(U[j], TmpPara);
			} else {
				out[i][j] = PutGamma(U[j], TmpPara);
			}
		}
		
	}
	return(out);
}

// Call-Put Partity Functions
double EUoption::PriceUsingPutCallParity (double U) const
{	
	OptionPara para = this->GetPara();
	string optiontype = this->GetType();
	if (optiontype == "C")
		return(CallPrice(U, para) + para.K*exp(-para.r*para.T) - U);
	else
		return (PutPrice(U, para) + U - para.K*exp(-para.r*para.T));
}

void EUoption::CheckPutCallParity(double U) const
{
	OptionPara para = this->GetPara();
	// float is employed here to get off the round error
	if ((float)(CallPrice(U, para) + para.K*exp(-para.r*para.T)) == (float)(PutPrice(U, para) + U)){
		cout << "The Call-Put Parity is hold." << endl;
	} else {
		cout << "The Call-Put Parity is not hold." << endl;
	}
}

//Exact Vs Numerical
vector<double> EUoption::ExactVSNumericalDelta(double U, const vector<double>& h) const
{
	int row_size = h.size();
	vector<double> out(row_size);
	double v1;
	double v2;
	double vE;
	for (int i = 0; i < row_size; i++)
	{	
		
			v2 = Price(U + h[i]);
			v1 = Price(U - h[i]);
			vE = Delta(U);
			out[i] = (v2-v1)/(2*h[i]) - vE;
	}
	return (out);
}


vector<double> EUoption::ExactVSNumericalGamma(double U, const vector<double>& h) const
{
	int row_size = h.size();
	vector<double>  out(row_size);
	double v1;
	double v2;
	double v3;
	double vE;
	for (int i = 0; i < row_size; i++)
	{	
		v1 = Price(U - h[i]);
		v2 = Price(U);
		v3 = Price(U + h[i]);
		vE = Gamma(U);
		out[i] = (v3-2*v2+v1)/(h[i]*h[i]) - vE;
		
	}
	return (out);
}

