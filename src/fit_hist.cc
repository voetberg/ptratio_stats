#include <iostream>           
#include <vector>             
#include <stdlib.h>           
#include <utility>            
#include <chrono>             
#include <complex> 

#include <TFile.h>            
#include <TChain.h>           
#include <TH1.h>              
#include <TMinuit.h>                               
                             

using namespace std;          

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

//Set Bins  
Int_t nbins = 38; 
vector<Double_t> hist_data; 
Int_t npar = 4; 

//Creates a function to fit to 
//3rd deg polynomial
Double_t funct(Double_t x, Double_t *par){
  Double_t val = par[0] + par[1]*x + par[2]*x*x + par[3]*x*x*x; 
  return val; 
}

//Returned the chi sq'd
void chisq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  //calculate chisquare
   Double_t chisq = 0;
   for (int i=0; i<nbins; i++) {
     Double_t calc = funct(hist_data[i],par); 
     chisq += pow((hist_data[i]-calc),2)/calc;
   }
   f=chisq;
}

int main(){ 
  const vector<string>
    num {"300","350","400","400p"}; 

  //Load in histogram 
  TFile *file = new TFile("/home/voetberg/voetberg/vari_dist/ptratio_data.root"); //Run from terminal
  //TFile *file = new TFile("$TMPDIR/ptratio_data.root"); //Run from Condor
  cout<<"Loaded in File"<<endl; 
  
  //Loop through each histogram in the file
  for (int i=0; i<4; ++i){
    TString name = "Data_Signal__ptratio_" + num[i]; 
    TH1D *h_fit = (TH1D*) file->Get(name); 
    cout<<"Loaded in histogram: "<<h_fit->GetName()<<endl;

    //Set the data in a container
    hist_data.reserve(nbins); 
    //const Double_t bin_width = 20./nbins;
    for (int i=1; i<=nbins; ++i){
      Double_t x = h_fit->GetBinContent(i); 
      hist_data.push_back(x); 
    }
    cout<<"Set data in vector"<<endl; 


    //Use Minuit
    TMinuit *gMinuit = new TMinuit(nbins); 
    gMinuit->SetFCN(chisq); 
    cout<<"Set Minuit"<<endl;  

    Double_t arglist[10]; 
    Int_t ier = 0; 

    arglist[0]=1; 
    gMinuit->mnexcm("SET ERR", arglist,1,ier); 
    
    //Set initial values 
    static Double_t param_start[4] ={100,-10,10,-1}; 
    static Double_t step = .01;
    gMinuit->mnparm(0, "a1", param_start[0], step, 0,0, ier); 
    gMinuit->mnparm(1, "a2", param_start[1], step, 0,0, ier); 
    gMinuit->mnparm(2, "a3", param_start[2], step, 0,0, ier); 
    gMinuit->mnparm(3, "a4", param_start[3], step, 0,0, ier); 
    cout<<"Set intial parameters"<<endl; 

    //Use migrad to minimize 
    arglist[0]=500; 
    arglist[1]=1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ier); 

    //Then Print
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnprin(3,amin);
     
  }

}
