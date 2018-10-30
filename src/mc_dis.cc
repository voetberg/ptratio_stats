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
                              
#include "Higgs2diphoton.hh"  
                              
using namespace std;          
bool sort_pt (TLorentzVector i, TLorentzVector j){return (i.Pt()>j.Pt());}                         
Higgs2diphoton Hdecay;        
                              
int main(){                                                               
  //===============================================                         
  TChain chain_s("t3"); 
  TChain chain_b("t3"); 

  chain_s.Add("$TMPDIR/H1.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_100.root",0); 
  chain_s.Add("$TMPDIR/H1.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_101.root",0); 
  
  chain_b.Add("$TMPDIR/born6.root",0); 
  chain_b.Add("$TMPDIR/born7.root",0); 
  chain_b.Add("$TMPDIR/born8.root",0); 
  chain_b.Add("$TMPDIR/born9.root",0); 
  
  cout<<"Read in files"<<endl; 
  
  //===============================================                            
  constexpr size_t max = 10;                                                    
  Float_t px_s[max], py_s[max], pz_s[max], E_s[max];                            
  Int_t nparticle_s, kf_s[max];                                                 
  Double_t weight_s;                                                            
                                                                                
  Float_t px_b[max], py_b[max], pz_b[max], E_b[max];                            
  Int_t nparticle_b, kf_b[max];                                                 
  Double_t weight_b;                                                            
  //===============================================
  chain_s.SetBranchAddress("nparticle",&nparticle_s);                          
  chain_s.SetBranchAddress("kf",kf_s);                                         
  chain_s.SetBranchAddress("px",px_s);                                         
  chain_s.SetBranchAddress("py",py_s);                                         
  chain_s.SetBranchAddress("pz",pz_s);                                         
  chain_s.SetBranchAddress("E",E_s);                                           
  chain_s.SetBranchAddress("weight2",&weight_s);                               
 
  chain_b.SetBranchAddress("nparticle",&nparticle_b);                          
  chain_b.SetBranchAddress("kf",kf_b);                                         
  chain_b.SetBranchAddress("px",px_b);                                         
  chain_b.SetBranchAddress("py",py_b);                                         
  chain_b.SetBranchAddress("pz",pz_b);                                         
  chain_b.SetBranchAddress("E",E_b);                                           
  chain_b.SetBranchAddress("weight2",&weight_b);                               
  //==============================================                              
  //Loading in data
  TChain chain_d("mini"); 

  chain_d.Add("$TMPDIR/data15.root"); 
  chain_d.Add("$TMPDIR/data16_DS1.root"); 
  chain_d.Add("$TMPDIR/data16_DS2.root"); 

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
  num.emplace_back("300");
  num.emplace_back("350");
  num.emplace_back("400");
  num.emplace_back("400p");

  for (int i=0; i<4; ++i){
    TString outname = "ptratio_dis_" + num[i] + ".root"; 
    TFile* out = TFile::Open(outname, "RECREATE"); 
    cout<<"Writing "<<outname<<"....."<<endl; 

    //===============================================                             
    //Ouputs
    int nbins = 38.;
    
    //MC
    //Signal
    TH1D* h_ptratio_s = new TH1D("Signal__ptratio_"+num[i], "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", nbins, 1, 20);                         //Sideband                
    TH1D* h_ptratio_b = new TH1D("Background__ptratio_"+num[i], "p_{T}^{#gamma_{1}} /p_{T}^{#gamma_{2}}", nbins, 1, 20);
                                                                
    //Data
    //Signal
    TH1D* h_ptratio_s_data = new TH1D("Data_Signal__ptratio_"+num[i], "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", nbins, 1, 20);
    //Sideband
    TH1D* h_ptratio_b_right = new TH1D("Data_Right__ptratio_"+num[i], "p_{T}^{#gamma_{1}} /p_{T}^{#gamma_{2}}", nbins, 1, 20);
    TH1D* h_ptratio_b_left = new TH1D("Data_Left__ptratio_"+num[i], "p_{T}^{#gamma_{1}} /p_{T}^{#gamma_{2}}", nbins, 1, 20);
    //==============================================                              
    const Long64_t entries_s = chain_s.GetEntries();                             
    //=============================================== 
    //Signal                                                                      
    for (long j=0; j<entries_s; ++j){                                             
      chain_s.GetEntry(j);                                                       
      TLorentzVector yy, y1, y2, jet;                                
      vector<TLorentzVector> jets;                                                
      for (long i=0; i<nparticle_s; ++i){                                         
        if (kf_s[i]==25){                                                         
          yy.SetPxPyPzE(px_s[i],py_s[i],pz_s[i],E_s[i]);                       
        }                                                                         
        else{                                                                     
          jets.emplace_back(px_s[i],py_s[i],pz_s[i],E_s[i]);                      
        }                                                                         
      }                                                                           
                                                                                  
      //Decay yy                                                               
      pair<TLorentzVector, TLorentzVector> diphoton = Hdecay(yy);              
                                                                                  
      //Sort Photons                                                                            
      if (diphoton.first.Pt()>diphoton.second.Pt()){                              
        y1 = diphoton.first;                                                 
        y2 = diphoton.second;                                                
      }                                                                           
      else{                                                                       
        y1 = diphoton.second;                                                
        y2 = diphoton.first;                                                 
      }                                                                           
                                                                                  
      jet=jets[0]; 
      
      //Cuts
      //Mass Cut                                                                  
      bool select = ((yy.M()>121.)&&(yy.M()<129.));                         
      //Rapidity Cut                                                              
      select &= (abs(y2.Rapidity())<2.4);                                    
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
      select &= (y2.DeltaR(jet)>.4);   
      //Jet Cuts                            
      select &= (jet.Pt()>30);              
      select &= (jet.Rapidity()<4.4);       
     
  
        //Diy Pt cuts 
      if (i==0){select &= (yy.Pt()>250. && yy.Pt()<300.);}                      
      if (i==1){select &= (yy.Pt()>300. && yy.Pt()<350);}                       
      if (i==2){select &= (yy.Pt()>350. && yy.Pt()<400);}                       
      if (i==3){select &= (yy.Pt()>400.);} 
      
      if (select){ 
        ptratio = (abs(y1.Pt())/abs(y2.Pt()));                          
        h_ptratio_s->Fill(ptratio, weight_s);                                           
      }       
    }
    //Background 
    //=============================================
    const Long64_t entries_b = chain_b.GetEntries(); 
    //=============================================
    for (long k=0; k<entries_b; ++k){
      chain_b.GetEntry(k); 

      vector<TLorentzVector> photons;                                           
      TLorentzVector y1, y2, jet;                                               
      photons.clear();                                                          
                                                                                
      //Assign particles                                                           
      for (long j=0; j<nparticle_b; ++j){                                       
        //if photon                                                                
        if (kf_b[j]== 22){                                                      
          photons.emplace_back(px_b[j], py_b[j], pz_b[j], E_b[j]);              
        }                                                                       
        //if jet                                                                   
        else{                                                                   
          jet.SetPxPyPzE(px_b[j], py_b[j], pz_b[j], E_b[j]);                    
        }                                                                       
      }                                                                         
      //Setting photons                                                            
      //Sort by pt                                                                 
      if (photons[0].Pt()>photons[1].Pt()){                                     
        y1 = photons[0];                                                        
        y2 = photons[1];                                                        
      }                                                                         
      else{                                                                     
        y1 = photons[1];                                                        
        y2 = photons[0];                                                        
      }                                                                         
      TLorentzVector yy = y1 + y2;  
      
      //Cuts
      //Mass Cut                                                                  
      bool select = ((yy.M()>121.)&&(yy.M()<129.));                         
      //Rapidity Cut                                                              
      select &= (abs(y2.Rapidity())<2.4);                                    
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
      select &= (y2.DeltaR(jet)>.4);   
      //Jet Cuts                            
      select &= (jet.Pt()>30);              
      select &= (jet.Rapidity()<4.4);       
     
  
        //Diy Pt cuts 
      if (i==0){select &= (yy.Pt()>250. && yy.Pt()<300.);}                      
      if (i==1){select &= (yy.Pt()>300. && yy.Pt()<350);}                       
      if (i==2){select &= (yy.Pt()>350. && yy.Pt()<400);}                       
      if (i==3){select &= (yy.Pt()>400.);} 

      
      if (select){ 
        ptratio = (abs(y1.Pt())/abs(y2.Pt()));                          
        h_ptratio_b->Fill(ptratio, weight_b);                                           
      }
    }

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
      //Jet Cuts                                                                
      for (int i=0; i<jet_n; ++i){                                              
        select &= (t_jet[i].Pt()>30);                                           
        select &= (t_jet[i].Rapidity()<4.4);                                    
      }              
      
      if (i==0){select &= (yy.Pt()>250. && yy.Pt()<300.);}                      
      if (i==1){select &= (yy.Pt()>300. && yy.Pt()<350);}                       
      if (i==2){select &= (yy.Pt()>350. && yy.Pt()<400);}                       
      if (i==3){select &= (yy.Pt()>400.);} 
      
      if (select){
        ptratio = abs(y1.Pt()/y2.Pt());

        bool sig = ((yy.M()>121.)&&(yy.M()<129.));
        bool left = ((yy.M()>105.)&&(yy.M()<121.));
        bool right = ((yy.M()>129.)&&(yy.M()<160.));
       
        if (sig){h_ptratio_s_data->Fill(ptratio);}
        if (left){h_ptratio_b_left->Fill(ptratio);}
        if (right){h_ptratio_b_right->Fill(ptratio);}
      }
    }

    cout<<outname<<" written"<<endl; 

    out->Write(); 
    out->Close(); 
  }
}
