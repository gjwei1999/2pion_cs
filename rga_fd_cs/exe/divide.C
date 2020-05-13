/****************************************/
/*                                      */
/*       Created by Jiawei Guo          */
/*           Duke University            */
/*                                      */
/****************************************/

#include "Lcore.h"

int main(const int argc, const char * argv[]){

    if (argc < 2){
    std::cout << "./divide <rebin ngroup> " << endl;
    return 0;
  }

    double density = 0.0708; //g/cm3
    double thickness = 5;//cm
    double mp = 1.6726219e-24;//g
    double e_charge = 1.60217662e-19; //C
    double convert;
    convert = 1.1494 * 0.95 * 1e-2 / e_charge * density / mp * thickness; //

    gStyle->SetPalette(53);
   
    int ngroup = atoi(argv[1]);
    float n_bin=500 / ngroup; 
    TString filename = "cs_rga1.root";
    
    TFile * fs = new TFile(filename.Data(), "RECREATE");
    
    //data
    TFile * fs_data = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/output3_0052_rebin.root", "r");
    
    TH1D * hdata_0 = (TH1D *) fs_data->Get("invmass_Ppip_twoPi");
    TH1D * hdata_1 = (TH1D *) fs_data->Get("invmass_Ppim_twoPi");
    TH1D * hdata_2 = (TH1D *) fs_data->Get("invmass_pippim_twoPi");
    TH1D * hdata_4 = (TH1D *) fs_data->Get("Q2_twoPi");
    TH2D * hdata_5 = (TH2D *) fs_data->Get("W_vs_q2_twoPi");
    TH2D * hdata_5b = (TH2D *) fs_data->Get("W_vs_q2_photon_twoPi");
    
    TH1D * hdata_6a = (TH1D *) fs_data->Get("W_q2_15_twoPi");
    TH1D * hdata_6b = (TH1D *) fs_data->Get("W_q2_2_twoPi");
    TH1D * hdata_6c = (TH1D *) fs_data->Get("W_q2_23_twoPi");
    TH1D * hdata_6d = (TH1D *) fs_data->Get("W_q2_34_twoPi");

    TH1D * hdata_7a = (TH1D *) fs_data->Get("W_q2_15_photon_twoPi");
    TH1D * hdata_7b = (TH1D *) fs_data->Get("W_q2_2_photon_twoPi");
    TH1D * hdata_7c = (TH1D *) fs_data->Get("W_q2_23_photon_twoPi");
    TH1D * hdata_7d = (TH1D *) fs_data->Get("W_q2_34_photon_twoPi");


    //simulation
    TFile * fs_simu = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/output3_simu_rebin.root", "r");
    TH1D * hsimu_0 = (TH1D *) fs_simu->Get("invmass_Ppip_twoPi");
    TH1D * hsimu_1 = (TH1D *) fs_simu->Get("invmass_Ppim_twoPi");
    TH1D * hsimu_2 = (TH1D *) fs_simu->Get("invmass_pippim_twoPi");
    TH1D * hsimu_4 = (TH1D *) fs_simu->Get("Q2_twoPi");
    TH2D * hsimu_5 = (TH2D *) fs_simu->Get("W_vs_q2_twoPi");
    TH2D * hsimu_5b = (TH2D *) fs_simu->Get("W_vs_q2_photon_twoPi");

    TH1D * hsimu_6a = (TH1D *) fs_simu->Get("W_q2_15_twoPi");
    TH1D * hsimu_6b = (TH1D *) fs_simu->Get("W_q2_2_twoPi");
    TH1D * hsimu_6c = (TH1D *) fs_simu->Get("W_q2_23_twoPi");
    TH1D * hsimu_6d = (TH1D *) fs_simu->Get("W_q2_34_twoPi");

    TH1D * hsimu_7a = (TH1D *) fs_simu->Get("W_q2_15_photon_twoPi");
    TH1D * hsimu_7b = (TH1D *) fs_simu->Get("W_q2_2_photon_twoPi");
    TH1D * hsimu_7c = (TH1D *) fs_simu->Get("W_q2_23_photon_twoPi");
    TH1D * hsimu_7d = (TH1D *) fs_simu->Get("W_q2_34_photon_twoPi");
    
    //original from generator
    TFile * fs_519 = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/original_574_rebin.root", "r");
    TH1D * h519_0 = (TH1D *) fs_519->Get("MpPip");
    TH1D * h519_1 = (TH1D *) fs_519->Get("MpPim");
    TH1D * h519_2 = (TH1D *) fs_519->Get("MPiPi");
    TH1D * h519_4 = (TH1D *) fs_519->Get("Q2");
    TH2D * h519_5 = (TH2D *) fs_519->Get("W_vs_q2");
    TH2D * h519_5b = (TH2D *) fs_519->Get("W_vs_q2_photon");
    
    TH1D * h519_6a = (TH1D *) fs_519->Get("W(q2<1.5)");
    TH1D * h519_6b = (TH1D *) fs_519->Get("W(q2<2)");
    TH1D * h519_6c = (TH1D *) fs_519->Get("W(2<q2<3)");
    TH1D * h519_6d = (TH1D *) fs_519->Get("W(3<q2<4)");
    TH2D * h519_7a = (TH2D *) fs_519->Get("W(q2<1.5)_photon");
    TH1D * h519_7b = (TH1D *) fs_519->Get("W(q2<2)_photon");
    TH1D * h519_7c = (TH1D *) fs_519->Get("W(2<q2<3)_photon");
    TH1D * h519_7d = (TH1D *) fs_519->Get("W(3<q2<4)_photon");
    
    TFile * fs_523 = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/original_575_rebin.root", "r");
    TH1D * h523_0 = (TH1D *) fs_523->Get("MpPip");
    TH1D * h523_1 = (TH1D *) fs_523->Get("MpPim");
    TH1D * h523_2 = (TH1D *) fs_523->Get("MPiPi");
    TH1D * h523_4 = (TH1D *) fs_523->Get("Q2");
    TH2D * h523_5 = (TH2D *) fs_523->Get("W_vs_q2");
    TH2D * h523_5b = (TH2D *) fs_523->Get("W_vs_q2_photon");

    TH1D * h523_6a = (TH1D *) fs_523->Get("W(q2<1.5)");
    TH1D * h523_6b = (TH1D *) fs_523->Get("W(q2<2)");
    TH1D * h523_6c = (TH1D *) fs_523->Get("W(2<q2<3)");
    TH1D * h523_6d = (TH1D *) fs_523->Get("W(3<q2<4)");
    TH2D * h523_7a = (TH2D *) fs_523->Get("W(q2<1.5)_photon");
    TH1D * h523_7b = (TH1D *) fs_523->Get("W(q2<2)_photon");
    TH1D * h523_7c = (TH1D *) fs_523->Get("W(2<q2<3)_photon");
    TH1D * h523_7d = (TH1D *) fs_523->Get("W(3<q2<4)_photon");
     
    TFile * fs_525 = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/original_576_rebin.root", "r");
    TH1D * h525_0 = (TH1D *) fs_525->Get("MpPip");
    TH1D * h525_1 = (TH1D *) fs_525->Get("MpPim");
    TH1D * h525_2 = (TH1D *) fs_525->Get("MPiPi");
    TH1D * h525_4 = (TH1D *) fs_525->Get("Q2");
    TH2D * h525_5 = (TH2D *) fs_525->Get("W_vs_q2");
    TH2D * h525_5b = (TH2D *) fs_525->Get("W_vs_q2_photon");

    TH1D * h525_6a = (TH1D *) fs_525->Get("W(q2<1.5)");
    TH1D * h525_6b = (TH1D *) fs_525->Get("W(q2<2)");
    TH1D * h525_6c = (TH1D *) fs_525->Get("W(2<q2<3)");
    TH1D * h525_6d = (TH1D *) fs_525->Get("W(3<q2<4)");
    TH2D * h525_7a = (TH2D *) fs_525->Get("W(q2<1.5)_photon");
    TH1D * h525_7b = (TH1D *) fs_525->Get("W(q2<2)_photon");
    TH1D * h525_7c = (TH1D *) fs_525->Get("W(2<q2<3)_photon");
    TH1D * h525_7d = (TH1D *) fs_525->Get("W(3<q2<4)_photon");

    TFile * fs_543 = new TFile("/work/clas12/users/gjwei/rga_fd_cs/after_rebin/original_577_rebin.root", "r");
    TH1D * h543_0 = (TH1D *) fs_543->Get("MpPip");
    TH1D * h543_1 = (TH1D *) fs_543->Get("MpPim");
    TH1D * h543_2 = (TH1D *) fs_543->Get("MPiPi");
    TH1D * h543_4 = (TH1D *) fs_543->Get("Q2");
    TH2D * h543_5 = (TH2D *) fs_543->Get("W_vs_q2");
    TH2D * h543_5b = (TH2D *) fs_543->Get("W_vs_q2_photon");

    TH1D * h543_6a = (TH1D *) fs_543->Get("W(q2<1.5)");
    TH1D * h543_6b = (TH1D *) fs_543->Get("W(q2<2)");
    TH1D * h543_6c = (TH1D *) fs_543->Get("W(2<q2<3)");
    TH1D * h543_6d = (TH1D *) fs_543->Get("W(3<q2<4)");
    TH2D * h543_7a = (TH2D *) fs_543->Get("W(q2<1.5)_photon");
    TH1D * h543_7b = (TH1D *) fs_543->Get("W(q2<2)_photon");
    TH1D * h543_7c = (TH1D *) fs_543->Get("W(2<q2<3)_photon");
    TH1D * h543_7d = (TH1D *) fs_543->Get("W(3<q2<4)_photon");

    //error histogram
    TH1D * herror_6c = new TH1D("W(2<q2<3)_error", "W(2<q2<3)_error", n_bin, 0, 5.5);
    TH1D * herror_7c = new TH1D("W(2<q2<3)_photon_error", "W(2<q2<3)_photon_error", n_bin, 0, 5.5);
    
    float error_6c = 0;
    float error_7c = 0;
    for(int i=1; i<n_bin + 1; i++){
        error_6c = hdata_6c->GetBinContent(i);
        error_7c = hdata_7c->GetBinContent(i);
        
        error_6c = sqrt(error_6c);
        error_7c = sqrt(error_7c);
        
        herror_6c->SetBinContent(i, error_6c);
        herror_7c->SetBinContent(i, error_7c);
    }

    h519_0->Add(h523_0,1);
    h519_0->Add(h525_0,1);
    h519_0->Add(h543_0,1);

    h519_1->Add(h523_1,1);
    h519_1->Add(h525_1,1);
    h519_1->Add(h543_1,1);
    
    h519_2->Add(h523_2,1);
    h519_2->Add(h525_2,1);
    h519_2->Add(h543_2,1);
    
    h519_4->Add(h523_4,1);
    h519_4->Add(h525_4,1);
    h519_4->Add(h543_4,1);
    
    h519_5->Add(h523_5,1);
    h519_5->Add(h525_5,1);
    h519_5->Add(h543_5,1);

    h519_5b->Add(h523_5b,1);
    h519_5b->Add(h525_5b,1);
    h519_5b->Add(h543_5b,1);

    h519_6a->Add(h523_6a,1);
    h519_6a->Add(h525_6a,1);
    h519_6a->Add(h543_6a,1);

    h519_6b->Add(h523_6b,1);
    h519_6b->Add(h525_6b,1);
    h519_6b->Add(h543_6b,1);

    h519_6c->Add(h523_6c,1);
    h519_6c->Add(h525_6c,1);
    h519_6c->Add(h543_6c,1);
    
    h519_6d->Add(h523_6d,1);
    h519_6d->Add(h525_6d,1);
    h519_6d->Add(h543_6d,1);
    
    h519_7a->Add(h523_7a,1);
    h519_7a->Add(h525_7a,1);
    h519_7a->Add(h543_7a,1);

    h519_7b->Add(h523_7b,1);
    h519_7b->Add(h525_7b,1);
    h519_7b->Add(h543_7b,1);

    h519_7c->Add(h523_7c,1);
    h519_7c->Add(h525_7c,1);
    h519_7c->Add(h543_7c,1);
    
    h519_7d->Add(h523_7d,1);
    h519_7d->Add(h525_7d,1);
    h519_7d->Add(h543_7d,1);
    
    //acceptance = simulation / original
    hsimu_0->Divide(h519_0);
    hsimu_1->Divide(h519_1);
    hsimu_2->Divide(h519_2);
    hsimu_4->Divide(h519_4);
    hsimu_5->Divide(h519_5);
    hsimu_5b->Divide(h519_5b);
    hsimu_6a->Divide(h519_6a);
    hsimu_6b->Divide(h519_6b);
    hsimu_6c->Divide(h519_6c);
    hsimu_6d->Divide(h519_6d);
    hsimu_7a->Divide(h519_7a);
    hsimu_7b->Divide(h519_7b);
    hsimu_7c->Divide(h519_7c);
    hsimu_7d->Divide(h519_7d);
    
    hsimu_0->SetDirectory(fs);
    hsimu_1->SetDirectory(fs);
    hsimu_2->SetDirectory(fs);
    hsimu_4->SetDirectory(fs);
    hsimu_5->SetDirectory(fs);
    hsimu_5b->SetDirectory(fs);
    hsimu_6a->SetDirectory(fs);
    hsimu_6b->SetDirectory(fs);
    hsimu_6c->SetDirectory(fs);
    hsimu_6d->SetDirectory(fs);
    hsimu_7a->SetDirectory(fs);
    hsimu_7b->SetDirectory(fs);
    hsimu_7c->SetDirectory(fs);
    hsimu_7d->SetDirectory(fs);
    
    hsimu_0->GetYaxis()->SetTitle("Acceptance");
    hsimu_1->GetYaxis()->SetTitle("Acceptance");
    hsimu_2->GetYaxis()->SetTitle("Acceptance");
    hsimu_4->GetYaxis()->SetTitle("Acceptance");
    hsimu_6a->GetYaxis()->SetTitle("Acceptance");
    hsimu_6b->GetYaxis()->SetTitle("Acceptance");
    hsimu_6c->GetYaxis()->SetTitle("Acceptance");
    hsimu_6d->GetYaxis()->SetTitle("Acceptance");
    hsimu_7a->GetYaxis()->SetTitle("Acceptance");
    hsimu_7b->GetYaxis()->SetTitle("Acceptance");
    hsimu_7c->GetYaxis()->SetTitle("Acceptance");
    hsimu_7d->GetYaxis()->SetTitle("Acceptance");
    hsimu_5->GetYaxis()->SetTitle("Q^2");
    hsimu_5b->GetYaxis()->SetTitle("Q^2");
    
    hsimu_4->GetXaxis()->SetTitle("W");
    hsimu_6a->GetXaxis()->SetTitle("W");
    hsimu_6b->GetXaxis()->SetTitle("W");
    hsimu_6c->GetXaxis()->SetTitle("W");
    hsimu_6d->GetXaxis()->SetTitle("W");
    hsimu_7a->GetXaxis()->SetTitle("W");
    hsimu_7b->GetXaxis()->SetTitle("W");
    hsimu_7c->GetXaxis()->SetTitle("W");
    hsimu_7d->GetXaxis()->SetTitle("W");

    hsimu_0->SetName("acceptance_invmass_Ppip");
    hsimu_1->SetName("acceptance_invmass_Ppim");
    hsimu_2->SetName("acceptance_invmass_pippim");
    hsimu_4->SetName("acceptance_Q2");
    hsimu_5->SetName("acceptance_WvsQ2");
    hsimu_5->SetName("acceptance_WvsQ2_photon");
    hsimu_6a->SetName("acceptance_elecprod_W(q2<1.5)");
    hsimu_6b->SetName("acceptance_elecprod_W(q2<2)");
    hsimu_6c->SetName("acceptance_elecprod_W(2<q2<3)");
    hsimu_6d->SetName("acceptance_elecprod_W(3<q2<4)");
    hsimu_7a->SetName("acceptance_photon_W(q2<1.5)");
    hsimu_7b->SetName("acceptance_photon_W(q2<2)");
    hsimu_7c->SetName("acceptance_photon_W(2<q2<3)");
    hsimu_7d->SetName("acceptance_photon_W(3<q2<4)");
//cross section
    hdata_0->Divide(hsimu_0);
    hdata_1->Divide(hsimu_1);
    hdata_2->Divide(hsimu_2);
    hdata_4->Divide(hsimu_4);
    hdata_5->Divide(hsimu_5);
    hdata_5b->Divide(hsimu_5b);
    hdata_6a->Divide(hsimu_6a);
    hdata_6b->Divide(hsimu_6b);
    hdata_6c->Divide(hsimu_6c);
    hdata_6d->Divide(hsimu_6d);
    hdata_7a->Divide(hsimu_7a);
    hdata_7b->Divide(hsimu_7b);
    hdata_7c->Divide(hsimu_7c);
    hdata_7d->Divide(hsimu_7d);

    herror_6c->Divide(hsimu_6c);
    herror_7c->Divide(hsimu_7c);

    hdata_0->Scale(1/convert);
    hdata_1->Scale(1/convert);
    hdata_2->Scale(1/convert);
    hdata_4->Scale(1/convert);
    hdata_5->Scale(1/convert);    
    hdata_5b->Scale(1/convert);
    hdata_6a->Scale(1/convert);
    hdata_6b->Scale(1/convert);
    hdata_6c->Scale(1/convert);
    hdata_6d->Scale(1/convert);
    hdata_7a->Scale(1/convert);
    hdata_7b->Scale(1/convert);
    hdata_7c->Scale(1/convert);
    hdata_7d->Scale(1/convert);

    herror_6c->Scale(1/convert);
    herror_7c->Scale(1/convert);

    //binwidth
    hdata_0->Scale(n_bin/(3.5-1.0));
    hdata_1->Scale(n_bin/(3.5-1.0));
    hdata_2->Scale(n_bin/2.5);
    hdata_4->Scale(n_bin/12.0);
    hdata_5->Scale(n_bin/12.0*n_bin/5.5);
    hdata_5b->Scale(n_bin/12.0*n_bin/5.5);
    hdata_6a->Scale(n_bin/5.5);
    hdata_6b->Scale(n_bin/5.5);
    hdata_6c->Scale(n_bin/5.5);
    hdata_6d->Scale(n_bin/5.5);
    hdata_7a->Scale(n_bin/5.5);
    hdata_7b->Scale(n_bin/5.5);
    hdata_7c->Scale(n_bin/5.5);
    hdata_7d->Scale(n_bin/5.5);

    herror_6c->Scale(n_bin/5.5);
    herror_7c->Scale(n_bin/5.5);
    
    //unit from cm2 to ub
    hdata_0->Scale(1e30);
    hdata_1->Scale(1e30);
    hdata_2->Scale(1e30);
    hdata_4->Scale(1e30); 
    hdata_5->Scale(1e30);
    hdata_5b->Scale(1e30);
    hdata_6a->Scale(1e30);
    hdata_6b->Scale(1e30);
    hdata_6c->Scale(1e30);
    hdata_6d->Scale(1e30);
    hdata_7a->Scale(1e30);
    hdata_7b->Scale(1e30);
    hdata_7c->Scale(1e30);
    hdata_7d->Scale(1e30);
    
    herror_6c->Scale(1e30);
    herror_7c->Scale(1e30);
    
    for(int i=1; i<n_bin + 1; i++){
        error_6c = herror_6c->GetBinContent(i);
        error_7c = herror_7c->GetBinContent(i);
        hdata_6c->SetBinError(i, error_6c);
        hdata_7c->SetBinError(i, error_7c);
    }

    hdata_0->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_1->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_2->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_4->GetYaxis()->SetTitle("#frac{d#sigma}{dQ^2} (#mu b / GeV^2)");
    hdata_5->GetYaxis()->SetTitle("Q^2 (GeV^2)");
    hdata_5->GetXaxis()->SetTitle("W (GeV)");
    hdata_5b->GetYaxis()->SetTitle("Q^2 (GeV^2)");
    hdata_5b->GetXaxis()->SetTitle("W (GeV)");
    hdata_6a->GetXaxis()->SetTitle("W (GeV)");
    hdata_6b->GetXaxis()->SetTitle("W (GeV)");
    hdata_6c->GetXaxis()->SetTitle("W (GeV)");
    hdata_6d->GetXaxis()->SetTitle("W (GeV)");
    hdata_7a->GetXaxis()->SetTitle("W (GeV)");
    hdata_7b->GetXaxis()->SetTitle("W (GeV)");
    hdata_7c->GetXaxis()->SetTitle("W (GeV)");
    hdata_7d->GetXaxis()->SetTitle("W (GeV)");
    
    hdata_6a->GetYaxis()->SetTitle("#frac{d#sigma}{dW}");
    hdata_6b->GetYaxis()->SetTitle("#frac{d#sigma}{dW}");
    hdata_6c->GetYaxis()->SetTitle("#frac{d#sigma}{dW}");
    hdata_6d->GetYaxis()->SetTitle("#frac{d#sigma}{dW}");
    hdata_7a->GetYaxis()->SetTitle("#sigma_{#nu}");
    hdata_7b->GetYaxis()->SetTitle("#sigma_{#nu}");
    hdata_7c->GetYaxis()->SetTitle("#sigma_{#nu}");
    hdata_7d->GetYaxis()->SetTitle("#sigma_{#nu}");
    
    hdata_5->SetName("WvsQ2");
    hdata_5b->SetName("WvsQ2_photon");
    hdata_6a->SetName("cs_elecprod_W(q2<1.5)");
    hdata_6b->SetName("cs_elecprod_W(q2<2)");
    hdata_6c->SetName("cs_elecprod_W(2<q2<3)");
    hdata_6d->SetName("cs_elecprod_W(3<q2<4)");
    hdata_7a->SetName("cs_photon_W(q2<1.5)");
    hdata_7b->SetName("cs_photon_W(q2<2)");
    hdata_7c->SetName("cs_photon_W(2<q2<3)");
    hdata_7d->SetName("cs_photon_W(3<q2<4)");

    hdata_7c->GetYaxis()->SetRangeUser(0.0, 3.0);

    //write
    hdata_0->SetDirectory(fs);
    hdata_1->SetDirectory(fs);
    hdata_2->SetDirectory(fs);
    hdata_4->SetDirectory(fs);
    hdata_5->SetDirectory(fs);    
    hdata_5b->SetDirectory(fs);
    hdata_6a->SetDirectory(fs);    
    hdata_6b->SetDirectory(fs);    
    hdata_6c->SetDirectory(fs);    
    hdata_6d->SetDirectory(fs);    
    hdata_7a->SetDirectory(fs);    
    hdata_7b->SetDirectory(fs);    
    hdata_7c->SetDirectory(fs);    
    hdata_7d->SetDirectory(fs);    
    fs->Write();
     
     return 1;
}
