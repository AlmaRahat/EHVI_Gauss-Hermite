//Command-line application for calculating the 3-Dimensional
//Expected Hypervolume Improvement.
//By Iris Hupkens, 2013
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
    if (p->f[0] >= P[i]->f[0] && p->f[1] >= P[i]->f[1]){
      cerr << "Individual " << (i+1) << " is dominated or the same as another point; removing." << endl;
      P.erase(P.begin()+i);
      nr++;
    }
  }
  return nr;
}

//Loads a testcase from the file with the name filename.
void loadtestcase(char *filename, deque<individual*> & testcase, double r[], double mu[], double s[]){
  ifstream file;
  int n, inds = 0;
  file.open(filename, ios::in);
  file >> n;
  for (int i=0;i<n;i++){
    individual * tempvidual = new individual;
    file >> tempvidual->f[0] >> tempvidual->f[1];
    tempvidual->f[0] = tempvidual->f[0];
    tempvidual->f[1] = tempvidual->f[1];
    checkdominance(testcase,tempvidual);
    testcase.push_back(tempvidual);
  }
  file >> r[0] >> r[1];
  file >> mu[0] >> mu[1];
  file >> s[0] >> s[1];
   file.close();
}

int main(int argc, char *argv[]){
  
  struct timeval before, after;
  
  int n;
  deque<individual*> testcase;
  double r[2];
  double mu[2];
  double s[2];
  
  vector<mus*> pdf;
  mus * tempmus = new mus;
  pdf.push_back(tempmus);
  cout << setprecision(10);
  if (argc > 1){
    loadtestcase(argv[1],testcase,r,mu, s);
    double answer = ehvi2d(testcase, r, mu, s);
    cout << answer << endl;
    return 0;
  }
}
