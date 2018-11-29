#include <iostream>           
#include <vector>             
#include <stdlib.h>           
#include <utility>            
#include <chrono>             
                              
#include <boost/optional.hpp> 
#include <TFile.h>            
#include <TChain.h>           
#include <TH1.h>              
#include <TLorentzVector.h>   
#include <TTree.h>            
#include <TTable.h> 

                              
using namespace std;          
bool sort_pt (TLorentzVector i, TLorentzVector j){return (i.Pt()>j.Pt());}                         
                              
int main(){                                                               
  //Loading in data
  TChain chain_d("mini"); 

  chain_d.Add("$TMPDIR/data1516.root"); 
  chain_d.Add("$TMPDIR/data17.root"); 

  //Setting chains 
  Int_t photon_n, jet_n;                                                        
  Float_t photon_m[max], photon_eta[max], photon_phi[max], photon_pt[max];      
  Float_t jet_m[max], jet_eta[max], jet_phi[max], jet_pt[max];                  
                                                                            
  chain_d.SetBranchAddress("photon_n", &photon_n);                                
  chain_d.SetBranchAddress("photon_m", photon_m);                                 
  chain_d.SetBranchAddress("photon_eta", photon_eta);                             
  chain_d.SetBranchAddress("photon_phi", photon_phi);                             
  chain_d.SetBranchAddress("photon_pt", photon_pt);                               
                                                                            
  chain_d.SetBranchAddress("jet_n", &jet_n);                                      
  chain_d.SetBranchAddress("jet_m", jet_m);                                       
  chain_d.SetBranchAddress("jet_eta", jet_eta);                                   
  chain_d.SetBranchAddress("jet_phi", jet_phi);                                   
  chain_d.SetBranchAddress("jet_pt", jet_pt);  

  //===============================================
  //ALL THE TREES                                                               
  Float_t ptratio;                                                              
  //===============================================                             
  vector<TString>num;
  num.emplace_back("350");
  num.emplace_back("400");
  num.emplace_back("400p");

  for (int i=0; i<3; ++i){
    TString outname = "ptratio_table_" + num[i] + ".root"; 
    TFile* out = TFile::Open(outname, "RECREATE"); 
    cout<<"Writing "<<outname<<"....."<<endl; 

    //===============================================                             
    //Ouputs
    int nbins = 38.;
    
    //Data
    //Signal
    TH1D* h_ptratio = new TH1D("Data_Signal__ptratio_"+num[i], "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", nbins, 1, 20);

    //Data 
    //Go through each event                                                     
    const Long64_t entries_d = chain_d.GetEntries();                                
    for (long k=0; k<entries_d; ++k){                                             
      chain_d.GetEntry(k);                                                        
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
      
      if (i==0){select &= (yy.Pt()>300. && yy.Pt()<350);}                       
      if (i==1){select &= (yy.Pt()>350. && yy.Pt()<400);}                       
      if (i==2){select &= (yy.Pt()>400.);} 
      
      if (select){
        ptratio = abs(y1.Pt()/y2.Pt());

        bool sig = ((yy.M()>121.)&&(yy.M()<129.));
        if (sig){h_ptratio->Fill(ptratio);}
      }
    }

    cout<<outname<<" written"<<endl; 

    //Print out the contents
    cout<<outname<<" signal contents"<<endl; 
    Float_t *bins = new Float[h_ptratio.GetSize()]; 
    for (int i=0; i<bins; ++i){
      cout<<i<<h_ptratio.GetBinContent(i)<<endl;      
    }


    out->Write(); 
    out->Close(); 
  } 
}
