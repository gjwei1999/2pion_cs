#include "Lcore.h"


int main(){
    //rga
    //492 493 498 499
    //rga_eft 
    //519 523 525 535(543)
    TString filename = "original_492.root";
    TFile * fs = new TFile(filename.Data(), "RECREATE");
  
    TLorentzVector ki[2], kf[5];
    TLorentzVector kprot;
    ki[0].SetXYZT(0, 0, sqrt(10.6041 * 10.6041 - PARTICLE::e.M() * PARTICLE::e.M()), 10.6041);
    kprot.SetXYZT(0, 0, 0, Mp);
    
    
    TH1D * h0a = new TH1D("MpPip", "MpPip", 500, 1.0, 3.5);
    TH1D * h0b = new TH1D("MpPim", "MpPim", 500, 1.0, 3.5);
    TH1D * h0c = new TH1D("MPiPi", "MPiPi", 500, 0.0, 2.5);
    TH1D * h0d = new TH1D("MissingSQ", "MissingSQ", 500, -0.2, 0.2);
    TH1D * h0e = new TH1D("Q2", "Q2", 500, 0.0, 12.0);
    
    TH2D * h0f = new TH2D("W_vs_q2","W_vs_q2",500,0.0,5.5,500,0.0,12.0);

    h0a->SetDirectory(fs);
    h0b->SetDirectory(fs);
    h0c->SetDirectory(fs);
    h0d->SetDirectory(fs);
    h0e->SetDirectory(fs);
    h0f->SetDirectory(fs);    
   
    gStyle->SetPalette(53);

    //path to lund file
    TString path1 = "/volatile/clas12/osg/gjwei/job_492/lund_dir/";       //path to lund files
    TString file = "out_rga";
  
    int file1 = 1;
    for(file1 ; file1 < 501; file1++){
      
    //TString inputfile = path1 + file + std::to_string(file1) + ".lund";
    string inputfile = Form("%s%s%d.lund", path1.Data(), file.Data(), file1);
    
    ifstream infile( inputfile.c_str() );
    
    double temp, px, py, pz, E, weight;
    TLorentzVector PP, Pa, Pb, Pc, Pd, Pe, W;
    while (infile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> weight){
        //k[0] electron   k[1]proton    k[2]pi+    k[3]pi-
        infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
        kf[0].SetXYZT(px, py, pz, E);
        infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
        kf[1].SetXYZT(px, py, pz, E);
        infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
        kf[2].SetXYZT(px, py, pz, E);
        infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
        kf[3].SetXYZT(px, py, pz, E);
        
        PP = kf[1] + kf[2] + kf[3];
        Pa = kf[1] + kf[2];
        Pb = kf[1] + kf[3];
        Pc = kf[2] + kf[3];
        Pd = ki[0] + kprot - kf[0] - PP;
        Pe = ki[0] - kf[0];
        W = ki[0] - kf[0] + kprot;
        h0a->Fill(Pa.M(), weight);
        h0b->Fill(Pb.M(), weight);
        h0c->Fill(Pc.M(), weight);
        h0d->Fill(Pd.M2(), weight);
        h0e->Fill(-Pe.M2(), weight);
        h0f->Fill(W.M(), -Pe.M2(), weight);
    }    
      infile.close();
      if(file1%100==0) {cout<<file1<<endl;}
  }

  fs->Write();
  return 1;
}
