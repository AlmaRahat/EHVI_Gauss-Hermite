#include "mex.h"   
#include <deque>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "ehvi_calculations.h"
#include "ehvi_sliceupdate.h"
#include "ehvi_montecarlo.h"
#include "string.h"
#include "ehvi_multi.h"
#include <vector>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
  
using namespace std;


 unsigned int
timeDif(
    struct timeval *before,
    struct timeval *after)
{
    return ((after->tv_sec * 1000000 + after->tv_usec) - 
        (before->tv_sec * 1000000 + before->tv_usec));
}


//Performs the EHVI calculations using the scheme requested by the user.
void doscheme(char *schemename, deque<individual*> & testcase, double r[], vector<mus*> & pdf){
  double answer;
  vector<double> answervector;
  if (pdf.size() == 1)
    if (strcmp(schemename,"sliceupdate") == 0){
      cerr << "Calculating with slice-update scheme..." << endl;
      answer = ehvi3d_sliceupdate(testcase, r, pdf[0]->mu, pdf[0]->s);
      cout << answer << endl;
    } else if (strcmp(schemename,"2term") == 0){
      cerr << "Calculating with 2-term scheme..." << endl;
      answer = ehvi3d_2term(testcase, r, pdf[0]->mu, pdf[0]->s);
      cout << answer << endl;
    } else if (strcmp(schemename,"5term") == 0){
      cerr << "Calculating with 5-term scheme..." << endl;
      answer = ehvi3d_5term(testcase, r, pdf[0]->mu, pdf[0]->s);
      cout << answer << endl;
    } else if (strcmp(schemename,"8term") == 0){
      cerr << "Calculating with 8-term scheme..." << endl;
      answer = ehvi3d_8term(testcase, r, pdf[0]->mu, pdf[0]->s);
      cout << answer << endl;
    } else if (strcmp(schemename,"montecarlo") == 0){
      cerr << "Calculating with Monte Carlo scheme (" << MONTECARLOS << " iterations)..." << endl;
      answer = ehvi3d_montecarlo(testcase, r, pdf[0]->mu, pdf[0]->s);
      cout << answer << endl;
    } else {
        cerr << "Scheme " << schemename << " does not exist. Proper options are:" << endl
             << "2term" << endl
             << "5term" << endl
             << "8term" << endl
             << "sliceupdate" << endl
             << "montecarlo" << endl;
    }
  else {
    if (strcmp(schemename,"sliceupdate") == 0){
      cerr << "Calculating with slice-update scheme (multi-version)..." << endl;
      answervector = ehvi3d_sliceupdate(testcase, r, pdf);
      for (int i=0;i<answervector.size();i++)
        cout << answervector[i] << endl;
    } else if (strcmp(schemename,"5term") == 0){
      cerr << "Calculating with 5-term scheme (multi-version)..." << endl;
      answervector = ehvi3d_5term(testcase, r, pdf);
      for (int i=0;i<answervector.size();i++)
        cout << answervector[i] << endl;
    } else {
        cerr << "Scheme " << schemename << " does not exist." << endl
             << "Multi-versions have only been implemented for the 5-term and slice-update schemes!" << endl;
    }
  }
}

//Checks if p dominates P. Removes points dominated by p from P and return the number of points removed.
int checkdominance(deque<individual*> & P, individual* p){
  int nr = 0;
  for (int i=P.size()-1;i>=0;i--){
    if (p->f[0] >= P[i]->f[0] && p->f[1] >= P[i]->f[1] && p->f[2] >= P[i]->f[2]){
      cerr << "Individual " << (i+1) << " is dominated or the same as another point; removing." << endl;
      P.erase(P.begin()+i);
      nr++;
    }
  }
  return nr;
}

//Loads a testcase from the file with the name filename.
void loadtestcase(char *filename, deque<individual*> & testcase, double r[], vector<mus*> & pdf){
  ifstream file;
  int n, inds = 0;
  file.open(filename, ios::in);
  file >> n;
  for (int i=0;i<n;i++){
    individual * tempvidual = new individual;
    file >> tempvidual->f[0] >> tempvidual->f[1] >> tempvidual->f[2];
	/* 
    tempvidual->f[0] = tempvidual->f[0];
    tempvidual->f[1] = tempvidual->f[1];
    tempvidual->f[2] = tempvidual->f[2]; */
	
    checkdominance(testcase,tempvidual);
    testcase.push_back(tempvidual);
  }
  file >> r[0] >> r[1] >> r[2];
  
  while (!file.eof()){
    if (inds > 0)
      pdf.push_back(new mus);
    file >> pdf[inds]->mu[0] >> pdf[inds]->mu[1] >> pdf[inds]->mu[2];
    file >> pdf[inds]->s[0] >> pdf[inds]->s[1] >> pdf[inds]->s[2];
	cout << pdf[inds]->mu[0] << " " << pdf[inds]->mu[1] << " " << pdf[inds]->mu[2] << " " << pdf[inds]->s[0] << " " << pdf[inds]->s[1] << " " << pdf[inds]->s[2] << endl;
    if (file.fail()){
      //We discover this while trying to read an individual and will end it here.
      pdf.pop_back();
      file.close();
      return;
    }
    inds++;
  }
  file.close();
}

vector<double>  loadtestcase_interface(int n,  double** data0, int m, double** data2, double* data1, deque<individual*> & testcase, double r[], vector<mus*> & pdf){

   int inds = 0;
   struct timeval before, after;
 
  for (int i=0;i<n;i++){
    individual * tempvidual = new individual;
    tempvidual->f[0] = *(*(data0 + i));
	tempvidual->f[1] = *(*(data0 + i) + 1);
	tempvidual->f[2] = *(*(data0 + i) + 2);
	
/* 	mexPrintf("The existing points are \n");
	mexPrintf("%f\n", tempvidual->f[0]);
	mexPrintf("%f\n", tempvidual->f[1]);
	mexPrintf("%f\n", tempvidual->f[2]); */
	
	
    checkdominance(testcase,tempvidual);
    testcase.push_back(tempvidual);
  }
  //file >> r[0] >> r[1] >> r[2];
  r[0] = *(data1);
  r[1] = *(data1 + 1);
  r[2] = *(data1 + 2);
  
  
/*    mexPrintf("Reference points\n");
   mexPrintf("%f\n", r[0]);
   mexPrintf("%f\n", r[1]);
   mexPrintf("%f\n", r[2]); */
  
  for (int i=0; i<m; i++){
    if (inds > 0)
		pdf.push_back(new mus);
	pdf[inds]->mu[0] = *(*(data2 + i));
	pdf[inds]->mu[1] = *(*(data2 + i) + 1);
	pdf[inds]->mu[2] = *(*(data2 + i) + 2);
	pdf[inds]->s[0] = *(*(data2 + i) + 3);
	pdf[inds]->s[1] = *(*(data2 + i) + 4);
	pdf[inds]->s[2] = *(*(data2 + i) + 5);
	
/* 	
	mexPrintf("Added points and sigma\n");
	mexPrintf("%f\n", pdf[inds]->mu[0]);
	mexPrintf("%f\n", pdf[inds]->mu[1]);
	mexPrintf("%f\n", pdf[inds]->mu[2]);
	mexPrintf("%f\n", pdf[inds]->s[0]);
	mexPrintf("%f\n", pdf[inds]->s[1]);
	mexPrintf("%f\n", pdf[inds]->s[2]); */
	
	
	
	
	//cout << pdf[inds]->mu[0] << " " << pdf[inds]->mu[1] << " " << pdf[inds]->mu[2] << " " << pdf[inds]->s[0] << " " << pdf[inds]->s[1] << " " << pdf[inds]->s[2] << endl;
   /*  if (file.fail()){
      //We discover this while trying to read an individual and will end it here.
      pdf.pop_back();
      file.close();
      return;
    } */
    inds++;
  }
  
   // gettimeofday(&before, NULL); 
	vector<double> ehvi_2d = ehvi3d_sliceupdate(testcase, r, pdf);
/* 	mexPrintf("The EHVI is: \n"); 
    for (int i=0;i<ehvi_2d.size();i++)
    mexPrintf("%f\n", ehvi_2d[i]);  */
	

	/* gettimeofday(&after, NULL);
	ofstream test("Time_IRS_update.txt", ios::app);
	test<<(timeDif(&before, &after))<<endl;  */
	
	return ehvi_2d;
}
  
  
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])  
{  
    #define p_exsit prhs[0]
    #define p_ref   prhs[1]
    #define p_eval  prhs[2]
    #define ehvi    plhs[0]
/* 	#define test_ref     plhs[1]
	#define test_eval    plhs[2] */
//  d_r_exi: dimension, the number of row
//  d_c_exi: dimension, the number of coloum 
    int d_r_exi,d_c_exi,d_r_eval,d_c_eval,d_r_r,d_c_r;
    double** data0, **data2; // data0: exist points; data2:evaluated points
    double*  data1;           // data1: reference points
    double *real_data_exi;
	double *real_data_ref;
	double *real_data_eval;
	double *temp;
	struct timeval before, after;
 
	deque<individual*> testcase;
	double r[DIMENSIONS];
	double mu[DIMENSIONS];
	double s[DIMENSIONS];
	vector<mus*> pdf;
	mus * tempmus = new mus;
	pdf.push_back(tempmus);
	cout << setprecision(10);
	

//	mexPrintf("Exisiting points\n");   
// Transfer real data from MATLAB to C++
// Transfer the existing points
    d_r_exi = mxGetM(p_exsit);
    d_c_exi = mxGetN(p_exsit);
    data0 = (double**)malloc(d_r_exi*sizeof(double*));
   
    real_data_exi = (double *)mxGetPr(p_exsit);
    for (int i=0;i<d_r_exi;i++){
		temp = real_data_exi;
        *(data0+i) = (double*)malloc(d_c_exi*sizeof(double));
        for (int j=0;j<d_c_exi;j++){
            (*(*(data0+i)+j)) = *(temp+j*d_r_exi);
            //mexPrintf("%f\n", (*(*(data0+i)+j)));
        }
		real_data_exi++;
    }
	
	//mexPrintf("Reference points\n"); 
// Transfer real data from MATLAB to C++
// Transfer the reference points
    d_r_r = mxGetM(p_ref);
    d_c_r = mxGetN(p_ref);
    data1 = (double*)malloc(d_c_r*sizeof(double));
	real_data_ref = (double *)mxGetPr(p_ref);
    for (int i=0;i<d_c_r;i++){
        *(data1 + i) = *real_data_ref++;
		//mexPrintf("%f\n", *(data1 + i)); 
    } 
	
//	mexPrintf("Evaluated points\n"); 
// Transfer real data from MATLAB to C++
// Transfer the evaluated points
	d_r_eval = mxGetM(p_eval);
    d_c_eval = mxGetN(p_eval);
    data2 = (double**)malloc(d_r_eval*sizeof(double*));
   
    real_data_eval = (double *)mxGetPr(p_eval);
    for (int i=0;i<d_r_eval;i++){
		temp = real_data_eval;
        *(data2+i) = (double*)malloc(d_c_eval*sizeof(double));
        for (int j=0;j<d_c_eval;j++){
            (*(*(data2+i)+j)) = *(temp+j*d_r_eval);
             //mexPrintf("%f\n", (*(*(data2+i)+j)));
        }
		real_data_eval++;
    }
	
	
	
	vector<double> answer = loadtestcase_interface(d_r_exi,data0, d_r_eval, data2, data1, testcase,r,pdf);
/* 	mexPrintf("Answer is \n"); 
	for (int i=0;i<answer.size();i++)
	mexPrintf("%f\n", answer[i]);  */
	
	
	
/* 	ehvi = mxCreateNumericMatrix(d_r_eval, 1, mxDOUBLE_CLASS, mxREAL);
    double* outputBuffer0 = (double*)mxGetData(ehvi); 
    for (int i=0;i<d_r_exi;i++){
            outputBuffer0[i] = answer[i];
    } */
    
	
	
	ehvi = mxCreateNumericMatrix(1, d_r_eval, mxDOUBLE_CLASS, mxREAL);
    double* outputBuffer1 = (double*)mxGetData(ehvi); 
    for (int i=0;i<d_r_eval;i++){
            outputBuffer1[i] = answer[i];
    }
	
	
	
	
// Output exist points
/*     ehvi = mxCreateNumericMatrix(d_r_exi, d_c_exi, mxDOUBLE_CLASS, mxREAL);
    double* outputBuffer0 = (double*)mxGetData(ehvi); 
    for (int i=0;i<d_r_exi;i++){
        for (int j=0;j<d_c_exi;j++){
            outputBuffer0[j + i*d_c_exi] = (*(*(data0+i)+j));
        }
    }
     */
	
	
// Output reference points	
/* 	test_ref = mxCreateNumericMatrix(d_r_r, d_c_r, mxDOUBLE_CLASS, mxREAL);
    double* outputBuffer1 = (double*)mxGetData(test_ref); 
    for (int i=0;i<d_r_r;i++){
        for (int j=0;j<d_c_r;j++){
            outputBuffer1[j + i*d_c_eval] = *(data1+i);
        }
    } */
	
	 
	 
// Output evaluated points
/* 	test_eval = mxCreateNumericMatrix(d_r_eval, d_c_eval, mxDOUBLE_CLASS, mxREAL);
    double* outputBuffer2 = (double*)mxGetData(test_eval); 
    for (int i=0;i<d_r_eval;i++){
        for (int j=0;j<d_c_eval;j++){
            outputBuffer2[j + i*d_c_eval] = (*(*(data2+i)+j));
        }
    }
     */
}  