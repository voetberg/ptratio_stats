#include <iostream>           
#include <vector>             
#include <stdlib.h>           
#include <utility>            
#include <chrono>             
#include <complex> 

#include <boost/optional.hpp> 
#include <TFile.h>            
#include <TChain.h>           
#include <TH1.h>              
#include <TLorentzVector.h>   
#include <TTree.h>            
#include <TMinuit.h>                               
                             

using namespace std;          
bool sort_pt (TLorentzVector i, TLorentzVector j){return (i.Pt()>j.Pt());}                         

/*
//Just a Legendre function 
//Yanked from https://github.com/ivankp/tt_phase/blob/63e64a3672e8fe364088566bfbb2474057fa3538/include/Legendre.hh
double Legendre(double x, const double* c) {
  const double x2 = x*x, x4 = x2*x2, x6 = x4*x2;

  const double p2 = 1.5*x2 - 0.5;
  const double p4 = 4.375*x4 - 3.75*x2 + 0.375;
  const double p6 = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125;

  const double c0 = sqrt(
    0.5 - (0.2*pow(c[0],2) + (1./9.)*pow(c[1],2) + (1./13.)*pow(c[2],2)) );

  const auto phase = polar<double>(1.,c[3]);

  return norm( c0 + c[0]*phase*p2 + c[1]*p4 + c[2]*p6 );
}
*/

//Create the containers for the values 





//Creates a function to fit to 
//3rd deg polynomial
Double_t funct(Double_t *x, Double_t *par){
  Double_t val = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0]; 
  return val; 
}

//Returned the chi sq'd
void chisq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   const Int_t nbins = 38; //I just know this 
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     delta  = (z[i]-funct(x[i],y[i],par))/errorz[i];
     chisq += delta*delta;
   }
   f = chisq;
}

int main(){ 
  const vector<string>
    num {300,350,400,400p}; 

  //Load in histogram 
  TFile *file = new TFile("$TMPDIR/ptratio_data.root"); 
  cout>>"Loaded in File">>endl; 
  
  //Loop through each histogram in the file
  for (int i=0, i<4, ++i){
    TH1D *h_fit = (TH1F*) file->Get("Data_Signal__ptratio_"+num[i]); 
    cout>>"Loaded in histogram: ">>h_fit->GetName()>>endl;
    
    //Use Minuit
    TMinuit *Minuit = new TMinuit(38); 
    gMinuit->SetFCN(chisq); 
    cout>>"Set Minuit">>endl;  

    Double_t arglist[10]; 
    Int_t ier = 0; 

    arglst[0]=1; 
    gMinuit->mnexcm("SET ERR", arglist,1,ier); 
    
    //Set initial values 
    static Double_t parm_start{50000,20000,20000,-1000}; 
    static Double_t step = .01;
    gMinuit->mnparm(0, "a1", param_start[0], step, 0,0, ier); 
    gMinuit->mnparm(1, "a2", param_start[1], step, 0,0, ier); 
    gMinuit->mnparm(2, "a3", param_start[2], step, 0,0, ier); 
    gMinuit->mnparm(3, "a4", param_start[3], step, 0,0, ier); 
    cout>>"Set intial parameters">>endl; 

    //Use migrad to minimize 
    arglist[0]=500; 
    arglist[1]=1;
    gMinuit->mnexxcm("MIGRAD",arglist,2,ier); 

    //Then Print
    cout>>h_fit->GetName()>>endl; 
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnprin(3,amin);
     
  }

}
