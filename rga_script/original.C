#include "Lcore.h"


int main(){
    //rga
    //492 493 498 499
    //rga_eft 
    //519 523 525 535(543)
    //rga_new
    //574 575 576 577
    TString filename = "original_576.root";
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
    TH1D * h0g = new TH1D("W", "W", 500, 0.0, 5.5);

    TH1D * h1a = new TH1D("W(q2<1.5)","W(q2<1.5)",500,0.0,5.5);
    TH1D * h1b = new TH1D("W(q2<2)","W(q2<2)",500,0.0,5.5);
    TH1D * h1c = new TH1D("W(2<q2<3)","W(2<q2<3)",500,0.0,5.5);
    TH1D * h1d = new TH1D("W(3<q2<4)","W(3<q2<4)",500,0.0,5.5);
    
    TH1D * h2a = new TH1D("W(q2<1.5)_photon","W(q2<1.5)_photon",500,0.0,5.5);
    TH1D * h2b = new TH1D("W(q2<2)_photon","W(q2<2)_photon",500,0.0,5.5);
    TH1D * h2c = new TH1D("W(2<q2<3)_photon","W(2<q2<3)_photon",500,0.0,5.5);
    TH1D * h2d = new TH1D("W(3<q2<4)_photon","W(3<q2<4)_photon",500,0.0,5.5);
    h0a->SetDirectory(fs);
    h0b->SetDirectory(fs);
    h0c->SetDirectory(fs);
    h0d->SetDirectory(fs);
    h0e->SetDirectory(fs);
    h0f->SetDirectory(fs);    
    h0g->SetDirectory(fs);
    h1a->SetDirectory(fs);
    h1b->SetDirectory(fs);
    h1c->SetDirectory(fs);
    h1d->SetDirectory(fs);
    h2a->SetDirectory(fs);
    h2b->SetDirectory(fs);
    h2c->SetDirectory(fs);
    h2d->SetDirectory(fs);

    h0f->SetOption("COLZ");
    gStyle->SetPalette(53);

    //path to lund file
    TString path1 = "/volatile/clas12/osg/gjwei/job_576/lund_dir/";       //path to lund files
    TString file = "out_rga";
 
    int file1 = 1;
    for(file1 ; file1 < 501; file1++){
      
    //TString inputfile = path1 + file + std::to_string(file1) + ".lund";
    string inputfile = Form("%s%s%d.lund", path1.Data(), file.Data(), file1);
    
    ifstream infile( inputfile.c_str() );
    
    double temp, px, py, pz, E, weight, q2, w_value, epsilont, flux;
    double nu, theta_e;
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
        h0g->Fill(W.M(), weight);

        q2 = -Pe.M2();
        w_value = W.M();
	theta_e = kf[0].Theta();
	nu = 10.6041 - kf[0].E();
	epsilont = 1/(1 + 2*(1+nu*nu/q2) * TMath::Tan(theta_e/2)*TMath::Tan(theta_e/2));
	//flux = 1/4/137/M_PI/10.6041/10.6041/Mp/Mp;// 
	flux = w_value*(w_value*w_value - Mp*Mp)/Mp/Mp/10.6041/10.6041/4.0/M_PI;
	flux = flux/q2 /(1-epsilont)/137.0;
	//std::cout<<"flux   "<<flux<<std::endl;
	if(q2 < 1.5){
		h1a->Fill(w_value, weight);
		h2a->Fill(w_value, weight/flux);
	}
	if(q2 < 2.0){
		h1b->Fill(w_value, weight);
		h2b->Fill(w_value, weight/flux);
	}
	if(q2 < 3.0 && q2 > 2.0){
		h1c->Fill(w_value, weight);
		h2c->Fill(w_value, weight/flux);
	}
	if(q2 < 4.0 && q2 > 3.0){
		h1d->Fill(w_value, weight);
		h2d->Fill(w_value, weight/flux);
	}
    }    
      infile.close();
      if(file1%100==0) {cout<<file1<<endl;}
  }

  fs->Write();
  return 1;
}
