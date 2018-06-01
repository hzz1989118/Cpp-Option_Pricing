// PAMoption.cpp
//
// Constructors and Member Functions of PAMoption Class


#include "PAMoption.hpp"
#include <cmath>
#include <iostream>


// Kernel Functions (Haug)
double PAMoption::CallPrice(double U, const OptionPara& para) const
{
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double sig2 =  _sig*_sig;
	double y1 = .5 - (_b/sig2) + sqrt((_b/sig2 - .5) * (_b/sig2 - .5) + (2*_r/sig2));
	return ((_K/(y1 - 1)) * pow( ((y1-1)/y1) * (U/_K) ,y1));

}

double PAMoption::PutPrice(double U, const OptionPara& para) const
{
	double _K = para.K;
	double _sig = para.sig;
	double _r = para.r;
	double _b = para.b;
	double sig2 =  _sig*_sig;
	double y2 = .5 - (_b/sig2) - sqrt((_b/sig2 - .5) * (_b/sig2 - .5) + (2*_r/sig2));
	return ((_K/(1 - y2)) * pow( ((y2-1)/y2) * (U/_K) ,y2));
}


/////////////////////////////////////////////////////////////////////////////////////


void PAMoption::copy( const PAMoption& o2)
{

	this->SetPara(o2.GetPara());
	this->SetType(o2.GetType());
	
}

PAMoption::PAMoption(){
	ExactOption();
}

PAMoption::PAMoption(const OptionPara& para, string optionType){
	//specific constructor
	this->SetPara(para);
	this->SetType(optionType);
	if (this->GetType() == "c")
		this->SetType("C");
}

PAMoption::PAMoption(const PAMoption& o2)
{ // Copy constructor
	copy(o2);
}


PAMoption::~PAMoption()
{

}


PAMoption& PAMoption::operator = (const PAMoption& option2)
{

	if (this == &option2) return *this;

	copy (option2);

	return *this;
}

// Functions that calculate option price and sensitivities
double PAMoption::Price(double U) const
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
vector<double> PAMoption::Price(const vector<double>& U) const
{	
	int size_v = U.size();
	vector<double> out(size_v);
	for (int i = 0; i <  size_v; i++){
		out[i] = Price(U[i]);
	}
	return(out);
}

// Matrix Pricer
vector<vector<double> > PAMoption::Price(const vector<double>& U, const vector<vector<double> >& ParaMatrix) const
{	

	int row_size = ParaMatrix.size(); // row size of output matrix
	int num_para = ParaMatrix[0].size(); // col num of input matrix
	int col_size = U.size(); // col size of output matrix

	if (num_para != 5) //Note: You can put any number in column T; Can't be ignored
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

