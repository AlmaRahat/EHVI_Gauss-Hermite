//This is the implementation of variants of the 5term and slice-update scheme in which
//a vector of Gaussian PDFs is used instead of just one.
#include <vector>
#include <deque>
#include <algorithm>
#include <math.h>
#include "ehvi_hvol.h"
#include "ehvi_multi.h"
#include <iostream>
#include "helper.h"
#include "ehvi_consts.h"
#include <math.h>



// These Parameters are used for 2D problem.
double psi1[1002][1001];
double psi2[1002][1001];
double gausscdf1[1002][1001];
double gausscdf2[1002][1001]; 
double sum; 
char flag1[1002][1001];
char flag2[1002][1001]; 
char fpsi1[1002][1001];
char fpsi2[1002][1001]; 
double exipsi1[1002][1001];
double exipsi2[1002][1001];
double ex1, ex2; 
double exipsi_i1[1002][1001]; double exipsi_i2[1002][1001]; double exipsi_i3[1002][1001]; double exipsi_i4[1002][1001];
double exipsi_ii1[1002][1001]; double exipsi_ii2[1002][1001]; double exipsi_ii3[1002][1001]; double exipsi_ii4[1002][1001]; 


using namespace std;


vector<double> ehvi3d_5term(deque<individual*> P, double r[], vector<mus *> & pdf){

//5-term 3-dimensional ehvi calculation scheme. Subtracts 4 quantities off a rectangular volume.
  vector<double> answer; //The eventual answer
  int n = P.size(); //Holds amount of points.
  double Sminus; //Correction term for the integral.
  deque<individual*> Py, Pz; //P sorted by y/z coordinate
  sort(P.begin(), P.end(), ycomparator);
  for (int i=0;i<P.size();i++){
    Py.push_back(P[i]);
  }
  sort(P.begin(), P.end(), zcomparator);
  for (unsigned int i=0;i<P.size();i++){
    Pz.push_back(P[i]);
  }
  sort(P.begin(), P.end(), xcomparator);
  for (int i=0;i<pdf.size();i++){
    answer.push_back(0);
  }
  for (int z=0;z<=n;z++){
    for (int y=0;y<=n;y++){
      for (int x=0;x<=n;x++){
        double v[DIMENSIONS]; //upper corner of hypervolume improvement box
        double cl[DIMENSIONS], cu[DIMENSIONS]; //Boundaries of grid cells
        cl[0] = (x == 0 ? r[0] : P[x-1]->f[0]);
        cl[1] = (y == 0 ? r[1] : Py[y-1]->f[1]);
        cl[2] = (z == 0 ? r[2] : Pz[z-1]->f[2]);
        cu[0] = (x == n ? INFINITY : P[x]->f[0]);
        cu[1] = (y == n ? INFINITY : Py[y]->f[1]);
        cu[2] = (z == n ? INFINITY : Pz[z]->f[2]);
        //We have to find v. This naive implementation is O(n) per iteration.
        v[0] = r[0];
        v[1] = r[1];
        v[2] = r[2];
        bool dominated = false;
        for (unsigned int i=0;i<P.size();i++){
          if (P[i]->f[0] >= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] >= cu[2]){
            dominated = true;
            break;
          }
          else if (P[i]->f[0] <= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] >= cu[2]){
            if (P[i]->f[0] > v[0])
              v[0] = P[i]->f[0];
          }
          else if (P[i]->f[0] >= cu[0] && P[i]->f[1] <= cu[1] && P[i]->f[2] >= cu[2]){
            if (P[i]->f[1] > v[1])
              v[1] = P[i]->f[1];
          }
          else if (P[i]->f[0] >= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] <= cu[2]){
            if (P[i]->f[2] > v[2])
              v[2] = P[i]->f[2];
          }
        }
        if (dominated)
          continue; //Cell's contribution is 0.
        Sminus = hvol3d(Pz, v, cl);
        double xslice = calculateslice(P, v, cl, 0);
        double yslice = calculateslice(Py, v, cl, 1);
        double zslice = calculateslice(Pz, v, cl, 2);
      //And then we integrate.
        for (int i=0;i<pdf.size();i++){
            double psi1 = exipsi(v[0],cl[0],pdf[i]->mu[0],pdf[i]->s[0]) - exipsi(v[0],cu[0],pdf[i]->mu[0],pdf[i]->s[0]);
            double psi2 = exipsi(v[1],cl[1],pdf[i]->mu[1],pdf[i]->s[1]) - exipsi(v[1],cu[1],pdf[i]->mu[1],pdf[i]->s[1]);
            double psi3 = exipsi(v[2],cl[2],pdf[i]->mu[2],pdf[i]->s[2]) - exipsi(v[2],cu[2],pdf[i]->mu[2],pdf[i]->s[2]);

            double gausscdf1 = gausscdf((cu[0]-pdf[i]->mu[0])/pdf[i]->s[0]) - gausscdf((cl[0]-pdf[i]->mu[0])/pdf[i]->s[0]);
            double gausscdf2 = gausscdf((cu[1]-pdf[i]->mu[1])/pdf[i]->s[1]) - gausscdf((cl[1]-pdf[i]->mu[1])/pdf[i]->s[1]);
            double gausscdf3 = gausscdf((cu[2]-pdf[i]->mu[2])/pdf[i]->s[2]) - gausscdf((cl[2]-pdf[i]->mu[2])/pdf[i]->s[2]);
            double sum = (psi1*psi2*psi3) - (Sminus*gausscdf1*gausscdf2*gausscdf3);
            //gausscdf represents chance of a point falling within the range [cl,cu);
            //psi = partial expected improvement
            //so psi - (gausscdf * (cl - v)) = p's expected distance from cl
            sum -= (xslice * gausscdf2 * gausscdf3 * (psi1 - (gausscdf1 * (cl[0]-v[0]))));
            sum -= (yslice * gausscdf1 * gausscdf3 * (psi2 - (gausscdf2 * (cl[1]-v[1]))));
            sum -= (zslice * gausscdf1 * gausscdf2 * (psi3 - (gausscdf3 * (cl[2]-v[2]))));
            if (sum > 0)
              answer[i] += sum;
        }
      }
    }
  }
  return answer;
}

vector<double> ehvi3d_sliceupdate(deque<individual*> P, double r[DIMENSIONS], vector<mus*> & pdf){
  vector<double> answer;
  int n = P.size(); //Holds amount of points.

  double mu1;
  double mu2;
  double s1;
  double s2;
  int jj;
	 
   
	

  thingy *Pstruct; //2D array with information about the shape of the dominated hypervolume
  deque<specialind*> Px, Py; //P sorted by x/y/z coordinate with extra information.
  double cellength[DIMENSIONS] = {0};
  while (answer.size() < pdf.size())
    answer.push_back(0);
		
			
  sort(P.begin(), P.end(), xcomparator);
  //double answer = 0; //The eventual answer
  int k = P.size(); //Holds amount of points.
  #ifdef NAIVE_DOMPOINTS
  deque<individual*> dompoints; //For the old-fashioned way.
  #endif
  double Sminus; //Correction term for the integral.
  int Sstart = k-1, Shorizontal = 0;
  //See thesis for explanation of how the O(1) iteration complexity
  //is reached. NOTE: i = y = f[1], j = x = f[0]
  int psize=pdf.size();
  for (int i=0; i<psize ; i++) answer[i]=0.0;
  for(int jj=0;jj<=psize;jj++)
  {
	for(int ii=0;ii<=k;ii++)
	{
		flag1[ii][jj]=1;
		flag2[ii][jj]=1;
		fpsi1[ii][jj]=1;
		fpsi2[ii][jj]=1;	
	}
  }
 
  for (int i=0;i<=k;i++){
    Sminus = 0;
    Shorizontal = Sstart;
    for (int j=k-i;j<=k;j++){
      double fmax[2]; //staircase width and height
      double cl1, cl2, cu1, cu2, temp; //Boundaries of grid cells
	  double cl1_previous,cl2_previous,cu1_previous,cu2_previous;
      if (j == k){
        fmax[1] = r[1];
        cu1 = INFINITY;
      }
      else {
        fmax[1] = P[j]->f[1];
        cu1 = P[j]->f[0];
      }
      if (i == k){
        fmax[0] = r[0];
        cu2 = INFINITY;
      }
      else {
        fmax[0] = P[k-i-1]->f[0];
        cu2 = P[k-i-1]->f[1];
      }
      cl1 = (j == 0 ? r[0] : P[j-1]->f[0]);
      cl2 = (i == 0 ? r[1] : P[k-i]->f[1]);
      //Cell boundaries have been decided. Determine Sminus.
      #ifdef NAIVE_DOMPOINTS
      dompoints.clear();
      for (int m = 0; m < k; m++){
        if (cl1 >= P[m]->f[0] && cl2 >= P[m]->f[1]){
          dompoints.push_back(P[m]);
        }
      }
      Sminus = calculateS(dompoints, fmax);
      #else
      if (Shorizontal > Sstart){
        Sminus += (P[Shorizontal]->f[0] - fmax[0]) * (P[Shorizontal]->f[1] - fmax[1]);
      }
      Shorizontal++;
      #endif
      //And then we integrate.
	  
	  for (int ii=0;ii<psize;ii++){
			 mu1=pdf[ii]->mu[0];
			 mu2=pdf[ii]->mu[1];
			 s1=pdf[ii]->s[0];
			 s2=pdf[ii]->s[1];
			 
			// Check exipsi exist or not 
			 if (fpsi1[j][ii])//Not exist
			 {
				 exipsi_i1[j][ii] = s1*gausspdf((cl1-mu1)/s1);
				 exipsi_i2[j][ii] = gausscdf((cl1-mu1)/s1);
				 exipsi_i3[j][ii] = s1*gausspdf((cu1-mu1)/s1);
				 exipsi_i4[j][ii] = gausscdf((cu1-mu1)/s1);
				 fpsi1[j][ii] = 0;
			 }
			 psi1[j][ii] = exipsi_i1[j][ii] + ((fmax[0]-mu1) * exipsi_i2[j][ii]) 
			             -(exipsi_i3[j][ii] + ((fmax[0]-mu1) * exipsi_i4[j][ii]));
			 //psi1 <0; no need to continue  
			 if (psi1[j][ii] < 0) break;
			 		 
			 
			 // Check
			 if (fpsi2[i][ii])
			 {
				exipsi_ii1[i][ii] = s2*gausspdf((cl2-mu2)/s2);
				exipsi_ii2[i][ii] = gausscdf((cl2-mu2)/s2);
				exipsi_ii3[i][ii] = s2*gausspdf((cu2-mu2)/s2);
				exipsi_ii4[i][ii] = gausscdf((cu2-mu2)/s2);
			        fpsi2[i][ii] = 0;
			 }
			 psi2[i][ii] = exipsi_ii1[i][ii] + ((fmax[1]-mu2)*exipsi_ii2[i][ii])
				     -(exipsi_ii3[i][ii] + ((fmax[1]-mu2)*exipsi_ii4[i][ii]));
			 
			 if (psi2[i][ii] < 0) break;
			

			 // Same trick for gausscdf
       	                 if (flag1[j][ii])
		            {gausscdf1[j][ii] = gausscdf((cu1-mu1)/s1) - gausscdf((cl1-mu1)/s1);flag1[j][ii]=0;}
			 if (flag2[i][ii])
		             {gausscdf2[i][ii] = gausscdf((cu2-mu2)/s2) - gausscdf((cl2-mu2)/s2);flag2[i][ii]=0;}  
		    
			 sum = (psi1[j][ii]*psi2[i][ii]) - (Sminus*gausscdf1[j][ii]*gausscdf2[i][ii]);
			 if (sum > 0) {answer[ii] = answer[ii]  + sum;}
	   }       
    }
    Sstart--;
  }
  return answer;
}
