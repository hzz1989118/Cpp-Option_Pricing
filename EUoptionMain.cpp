// This Main cpp file provides all answers to questions in Group A and B.
// Comments along with code will indicate the questions

#include "EUoption.hpp"
#include "PAMoption.hpp"
#include <cmath>


int main()
{ 	


	////////////////////////////////////// Group A //////////////////////////////////////
	
	/////////////   Part1-Question(a)    /////////////////
	
	OptionPara EUopt1_para = {1.5, 120, .4, .04, 0};
	EUoption EUopt1(EUopt1_para, "C");
	cout << EUopt1.Price(100.0) << " ";
	EUopt1.toggle();
	cout << EUopt1.Price(100.0) << endl;
	EUopt1.toggle();

	cout << EUopt1.Delta(100.0) << endl;
	cout << EUopt1.Gamma(100.0) << endl;

	// OptionPara EUopt2_para = {1.0, 120, .4, .04, .0};
	// EUoption EUopt2(EUopt2_para, "C");
	// cout << EUopt2.Price(100) << " ";
	// EUopt2.toggle();
	// cout << EUopt2.Price(100) << endl;
	// EUopt2.toggle();

	// OptionPara EUopt3_para = {1.0, 10, .5, .12, .12};
	// EUoption EUopt3(EUopt3_para, "C");
	// cout << EUopt3.Price(5) << " ";
	// EUopt3.toggle();
	// cout << EUopt3.Price(5) << endl;
	// EUopt3.toggle();

	// OptionPara EUopt4_para = {30, 100, .3, .08, .08};
	// EUoption EUopt4(EUopt4_para, "C");
	// cout << EUopt4.Price(100) << " ";
	// EUopt4.toggle();
	// cout << EUopt4.Price(100) << endl;
	// EUopt4.toggle();

	// // //////////////////  Part1-Question(b)   /////////////////
	// cout << "For Batch 1, the put calculated by the Put-Call parity is " << EUopt1.PriceUsingPutCallParity(60.0) << endl;
	// EUopt1.CheckPutCallParity(60.0);

	// cout << "For Batch 2, the put calculated by the Put-Call parity is " << EUopt2.PriceUsingPutCallParity(100.0) << endl;
	// EUopt2.CheckPutCallParity(100.0);

	// cout << "For Batch 3, the put calculated by the Put-Call parity is " << EUopt3.PriceUsingPutCallParity(5.0) << endl;
	// EUopt3.CheckPutCallParity(5.0);

	// cout << "For Batch 4, the put calculated by the Put-Call parity is " << EUopt4.PriceUsingPutCallParity(100.0) << endl;
	// EUopt4.CheckPutCallParity(100.0);

	// // //////////////    Part1-Question(c)      ////////////
	// //initialize the vector
	// vector<double> vecU = getMesh(10, 50, 1.0);
	// cout << "\n\nS in vector: " << endl;
	// Print1d(vecU);
	
	// //calculate the vector price and print out
	// vector<double> C_vec_batch1 = EUopt1.Price(vecU);
	// //print
	// cout << "\nC of batch 1 in vector: " << endl;
	// Print1d(C_vec_batch1);

	// // //////////    Part1-Question(d)    /////////////

	// //initialize the matrix of S
	// vector<vector<double> > matP(2, vector<double> (5));
	// double arr1[] = {.25, 65.0, .30, .08, .08};
	// double arr2[] = {1.0, 100, .2, .0, .0};
	// for (int i; i < 5; i++){
	// 	matP[0][i] = arr1[i];
	// 	matP[1][i] = arr2[i];
	// }
	
	// cout << "\n\nInput in matrix: " << endl;
	// cout << matP[0].size() << endl;
	// Print2d(matP);

	// // Obtain matrix of C by matrix of S
	// vector<vector<double> > C_mat = OptionMatrixPricer(EUopt1, vecU , matP);
	// // // Print the matrix out
	// cout << "C of batch 1 in matrix : " << endl;
	// Print2d(C_mat);

	// // //////////    Part2-Question(a)    ////////////
	// // initialize option for batch 5
	// OptionPara EUopt5_para = {.5, 100, .36, .1, 0};
	// EUoption EUopt5(EUopt5_para, "C");
	// cout << "\nCall: " << EUopt5.Price(105);
	// EUopt5.toggle();
	// cout << " Put: " << EUopt5.Price(105) << endl;
	// EUopt5.toggle();

	// cout << "Gamma: " << EUopt5.Gamma(105.0) << endl;
	// cout << "Delta for Call: " << EUopt5.Delta(105) << " and ";
	// EUopt5.toggle(); // toggle back to call
	// cout << "for Put: " << EUopt5.Delta(105) << endl;
	// EUopt5.toggle(); // toggle back to call

	// // //////////    Part2-Question(b)    ////////////
	// // Delta for Batch 5 at S for 10 - 50;
	// vector<double> Delta_vec_Euopt5 = EUopt5.Price(vecU);
	// //print
	// cout << "\nDelta in vector: " << endl;
	// Print1d(Delta_vec_Euopt5);

	// // //////////    Part2-Question(c)    ////////////
	// // Gamma/Delta Matrix 
	// vector<vector<double> > Delta_mat = OptionMatrixDelta(EUopt5, vecU , matP);
	// // Print the matrix out
	// cout << "\n\nDelta in matrix : " << endl;
	// Print2d(Delta_mat); 

	// vector<vector<double> > Gamma_mat = OptionMatrixGamma(EUopt5, vecU , matP);
	// // Print the matrix out
	// cout << "Gamma in matrix : " << endl;
	// Print2d(Gamma_mat); 

	// // // //////////    Part2-Question(d)    ////////////
	// // //initialize the vector for possible S

	// vector<double> vech = getMesh(0.05, 5, 0.05);
	// cout << "\nh in vector: " << endl;
	// Print1d(vech);

	// vector<double> EvsN_Delta = EUopt5.ExactVSNumericalDelta(105, vech);
	// cout << "\nNumerical Delta - Exact Delta: " << endl;
	// Print1d(EvsN_Delta);

	// cout << "\nNumerical Gamma - Exact Gamma: " << endl;
	// vector<double> EvsN_Gamma = EUopt5.ExactVSNumericalGamma(105, vech);
	// Print1d(EvsN_Gamma);


	// ////////////////////////////////////// Group B //////////////////////////////////////
	// //////////    Question(a&b)    ////////////
	// OptionPara PAMopt1_para = {999.0, 100.0, .1, .1, 0.02};
	// PAMoption PAMopt1(PAMopt1_para, "C");
	// cout << "\n\nPerpetual American Option in Group B is with Call: " << PAMopt1.Price(110) << " ";
	// PAMopt1.toggle();
	// cout << "and Put: " << PAMopt1.Price(110) << endl;
	// PAMopt1.toggle();

	// //////////    Question(c&d)    ////////////
	// cout << "\n\nFor S in vector: " << endl; // S vector
	// Print1d(vecU);

	// vector<double> C_vec_Qc = PAMopt1.Price(vecU); //Call price vector
	// //print
	// cout << "\nthe corresponding Call: " << endl;
	// Print1d(C_vec_Qc);


	// vector<vector<double> > matPara(2, vector<double> (5));
	// double para1[] ={999.0, 100.0, .1, .1, 0.02};
	// double para2[] ={999.0, 60.0, .05, .08, 0.02};
	// for(int j ; j<5; j++){
	// 	matPara[0][j] = para1[j];
	// 	matPara[1][j] = para2[j];
	// }
	// cout << "\n\nInput parameters in matrix: " << endl; //Para Matrix
	// Print2d(matPara);

	// vector<vector<double> > C_mat_Qd = PAMopt1.Price(vecU , matPara); // Call Price Matrix
	// // Print the matrix out
	// cout << "\nThe corresponding call in matrix : " << endl;
	// Print2d(C_mat_Qd);


	return 0;
}