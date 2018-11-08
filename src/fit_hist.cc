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


auto hist(TString num){
  //Returns a histogram 
  //Takes num argument to decide pT cut range (~line 130)    
  //===============================================                            
  //Loading in data
  TChain chain("mini"); 

  chain.Add("$TMPDIR/data1516.root"); 
  chain.Add("$TMPDIR/data17.root"); 

  //Setting chains 
  constexpr size_t max = 10;                                                    
  Int_t photon_n, jet_n;                                                        
  Float_t photon_m[max], photon_eta[max], photon_phi[max], photon_pt[max];      
  Float_t jet_m[max], jet_eta[max], jet_phi[max], jet_pt[max];                  
                                                                            
  chain.SetBranchAddress("photon_n", &photon_n);                                
  chain.SetBranchAddress("photon_m", photon_m);                                 
  chain.SetBranchAddress("photon_eta", photon_eta);                             
  chain.SetBranchAddress("photon_phi", photon_phi);                             
  chain.SetBranchAddress("photon_pt", photon_pt);                               
                                                                            
  chain.SetBranchAddress("jet_n", &jet_n);                                      
  chain.SetBranchAddress("jet_m", jet_m);                                       
  chain.SetBranchAddress("jet_eta", jet_eta);                                   
  chain.SetBranchAddress("jet_phi", jet_phi);                                   
  chain.SetBranchAddress("jet_pt", jet_pt);  

  //===============================================
  //ALL THE TREES                                                               
  Float_t ptratio;                                                              
  //===============================================                             
  //Ouputs
  int nbins = 38.;
  TH1D* h_ptratio_s_data = new TH1D("Data_Signal__ptratio_"+num, "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", nbins, 1, 20);
  //==============================================                              
  const Long64_t entries = chain.GetEntries();                                
  for (long k=0; k<entries; ++k){                                             
    chain.GetEntry(k);                                                        
    //Assign Particles                                                        
    TLorentzVector yy, y1, y2, jet, photons, jets;                            
    vector<TLorentzVector> t_y, t_jet;                                        
                                                                             
    for (int j=0; j<photon_n; ++j){                                           
      photons.SetPtEtaPhiM(photon_pt[j], photon_eta[j], photon_phi[j], photon_m[j]); 
      t_y.emplace_back(photons);                                              
    }                                                                         
    if (jet_n!=0){                                                            
      for (int j=0; j<jet_n; ++j){                                            
        jets.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], jet_m[j]);       
        t_jet.emplace_back(jets);                                             
      }                                                                       
    }                                                                         
    else{                                                                     
      jets.SetPtEtaPhiM(0,0,0,0);                                             
      t_jet.emplace_back(jets);                                               
    }                                                                         
                                                                             
    //Sort by pT                                                              
    sort(t_jet.begin(), t_jet.end(), sort_pt);                                
    sort(t_y.begin(), t_y.end(), sort_pt);                                    
                                                                             
    jet = t_jet[0];                                                           
    y1 = t_y[0];                                                              
    y2 = t_y[1];                                                              
    yy = y1+y2;
    
    //Cuts                                                                    
    //Rapidity Cut                                                               
    bool select = (abs(y2.Rapidity())<2.4);                                   
    select &= (abs(y1.Rapidity())<2.4);                                       
    //PseudoRapidity Cut                                                         
    select &= (abs(y1.PseudoRapidity())<2.37);                                
    select &= !(1.37<abs(y1.PseudoRapidity()) && abs(y1.PseudoRapidity())<1.52);
    select &= (abs(y2.PseudoRapidity())<2.37);                                
    select &= !(1.37<abs(y2.PseudoRapidity()) && abs(y2.PseudoRapidity())<1.52);
    //Pt Cut                                                                     
    select &= (y1.Pt()>.35*yy.M());                                           
    select &= (y2.Pt()>.25*yy.M());                                           
    //Delta R                                                                    
    select &= (y1.DeltaR(jet)>.4);                                            
    select &= (y2.DeltaR(jet)>.4);                                            
    //Jet Cuts                                                                
    for (int i=0; i<jet_n; ++i){                                              
      select &= (t_jet[i].Pt()>30);                                           
      select &= (t_jet[i].Rapidity()<4.4);                                    
    }              
    
    if (num == "350"){select &= (yy.Pt()>300. && yy.Pt()<350);}                       
    if (num == "400"){select &= (yy.Pt()>350. && yy.Pt()<400);}                       
    if (num == "400p"){select &= (yy.Pt()>400.);} 
    
    if (select){
      ptratio = abs(y1.Pt()/y2.Pt());
      bool sig = ((yy.M()>122.)&&(yy.M()<128.));
      if (sig){h_ptratio_s_data->Fill(ptratio);}
    }
  }
  return (h_ptratio_s_data); 
}

//chi^2 function for fit
//I stole these from Ivan so ?????

//Just a Legendre function 
//Yanked from https://github.com/ivankp/tt_phase/blob/63e64a3672e8fe364088566bfbb2474057fa3538/include/Legendre.hh
//Used for the chi2'd
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


//Returns the chi2'd for a given histogram 
auto chi2_maker(int nbins,auto hist){
  vector<array<double,3>> chi2_data; 
  chi2_data.reserve(nbins); 
  const double bin_w = (hist.axis().max() - hist.axis().min())/nbins;
  
  for (unsigned i=0; i<nbins; ++i){
    const auto& b = hist.bins()[i]; 
    chi2_data.push_back({
        b.w/(bin_w), 
        b.w2/sqrt(bin_w),
        hist.axis().min() + (i+0.5)*bin_w
        }); 
  }
  double chi2 = 0; 
  for (unsigned i=0; i<nbins; ++i){
    chi2+= pow(b[i][0] - Legendre(b[i][2],c),2)/b[i][1]; 
  }
  return chi2; 
}

//Uses the chi'd and minuit to produce a fit
auto Chi2_Fit(){
  
}


int main(){ 
  vector<TString>num;
  num.emplace_back("350");
  num.emplace_back("400");
  num.emplace_back("400p");

  for(int i=0; i<3; ++i){
    //Get the histogram 
    auto histogram = hist(num[i]); 
    //-2 to account for underflow and overflow bins 
    int nbins = histogram->GetSize()-2;

     
  
  }
}
