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
    convert = 802692 * 0.95 * 1e-9 / e_charge * density / mp * thickness; //1.14e39  

    gStyle->SetPalette(53);
   
    int ngroup = atoi(argv[1]); 
    float n_bin = 500 / ngroup;
    
    TString filename = "cs_rga_eft1.root";
    
    TFile * fs = new TFile(filename.Data(), "RECREATE");
    
    //data
    TFile * fs_data = new TFile("/work/clas12/users/gjwei/FT_analysis/2pion_ftpid/build_new1/output3_ftpid_strict_cut.root", "r");
    
    TH1D * hdata_0 = (TH1D *) fs_data->Get("invmass_Ppip_twoPi");
    TH1D * hdata_1 = (TH1D *) fs_data->Get("invmass_Ppim_twoPi");
    TH1D * hdata_2 = (TH1D *) fs_data->Get("invmass_pippim_twoPi");
    TH1D * hdata_4 = (TH1D *) fs_data->Get("Q2_twoPi");
    
    TH2D * hdata_5 = (TH2D *) fs_data->Get("W_vs_q2_twoPi");
    TH2D * hdata_5b = (TH2D *) fs_data->Get("W_vs_q2_photon_twoPi");    

    TH1D * hdata_6a = (TH1D *) fs_data->Get("W_photon_twoPi");
    TH1D * hdata_6b = (TH1D *) fs_data->Get("W_twoPi");   

    TH1D * hdata_7a = (TH1D *) fs_data->Get("W_03_q2_05_twoPi");
    TH1D * hdata_7b = (TH1D *) fs_data->Get("W_03_q2_05_photon_twoPi");
     
    //simulation
    TFile * fs_simu = new TFile("/work/clas12/users/gjwei/FT_analysis/2pion_ftpid/build_new1/output3_ftpid_strict_cut_simu.root", "r");
    TH1D * hsimu_0 = (TH1D *) fs_simu->Get("invmass_Ppip_twoPi");
    TH1D * hsimu_1 = (TH1D *) fs_simu->Get("invmass_Ppim_twoPi");
    TH1D * hsimu_2 = (TH1D *) fs_simu->Get("invmass_pippim_twoPi");
    TH1D * hsimu_4 = (TH1D *) fs_simu->Get("Q2_twoPi");

    TH2D * hsimu_5 = (TH2D *) fs_simu->Get("W_vs_q2_twoPi");
    TH2D * hsimu_5b = (TH2D *) fs_simu->Get("W_vs_q2_photon_twoPi"); 

    TH1D * hsimu_6a = (TH1D *) fs_simu->Get("W_photon_twoPi");
    TH1D * hsimu_6b = (TH1D *) fs_simu->Get("W_twoPi");

    TH1D * hsimu_7a = (TH1D *) fs_data->Get("W_03_q2_05_twoPi");
    TH1D * hsimu_7b = (TH1D *) fs_data->Get("W_03_q2_05_photon_twoPi");
    
    //original from generator
    TFile * fs_519 = new TFile("/work/clas12/users/gjwei/rga_ft_cs/ori_root/original_519.root", "r");
    TH1D * h519_0 = (TH1D *) fs_519->Get("MpPip");
    TH1D * h519_1 = (TH1D *) fs_519->Get("MpPim");
    TH1D * h519_2 = (TH1D *) fs_519->Get("MPiPi");
    TH1D * h519_4 = (TH1D *) fs_519->Get("Q2");

    TH2D * h519_5 = (TH2D *) fs_519->Get("W_vs_q2");
    TH2D * h519_5b = (TH2D *) fs_519->Get("W_vs_q2_photon");

    TH1D * h519_6a = (TH1D *) fs_519->Get("W_photon");
    TH1D * h519_6b = (TH1D *) fs_519->Get("W");

    TH1D * h519_7a = (TH1D *) fs_519->Get("W(03<q2<05)");
    TH1D * h519_7b = (TH1D *) fs_519->Get("W(03<q2<05)_photon");

    TFile * fs_523 = new TFile("/work/clas12/users/gjwei/rga_ft_cs/ori_root/original_523.root", "r");
    TH1D * h523_0 = (TH1D *) fs_523->Get("MpPip");
    TH1D * h523_1 = (TH1D *) fs_523->Get("MpPim");
    TH1D * h523_2 = (TH1D *) fs_523->Get("MPiPi");
    TH1D * h523_4 = (TH1D *) fs_523->Get("Q2");
    
    TH2D * h523_5 = (TH2D *) fs_523->Get("W_vs_q2");
    TH2D * h523_5b = (TH2D *) fs_523->Get("W_vs_q2_photon");  
  
    TH1D * h523_6a = (TH1D *) fs_523->Get("W_photon");
    TH1D * h523_6b = (TH1D *) fs_523->Get("W");
    
    TH1D * h523_7a = (TH1D *) fs_523->Get("W(03<q2<05)");
    TH1D * h523_7b = (TH1D *) fs_523->Get("W(03<q2<05)_photon");


    TFile * fs_525 = new TFile("/work/clas12/users/gjwei/rga_ft_cs/ori_root/original_525.root", "r");
    TH1D * h525_0 = (TH1D *) fs_525->Get("MpPip");
    TH1D * h525_1 = (TH1D *) fs_525->Get("MpPim");
    TH1D * h525_2 = (TH1D *) fs_525->Get("MPiPi");
    TH1D * h525_4 = (TH1D *) fs_525->Get("Q2");

    TH2D * h525_5 = (TH2D *) fs_525->Get("W_vs_q2");
    TH2D * h525_5b = (TH2D *) fs_525->Get("W_vs_q2_photon");

    TH1D * h525_6a = (TH1D *) fs_525->Get("W_photon");
    TH1D * h525_6b = (TH1D *) fs_525->Get("W");

    TH1D * h525_7a = (TH1D *) fs_525->Get("W(03<q2<05)");
    TH1D * h525_7b = (TH1D *) fs_525->Get("W(03<q2<05)_photon");


    TFile * fs_543 = new TFile("/work/clas12/users/gjwei/rga_ft_cs/ori_root/original_543.root", "r");
    TH1D * h543_0 = (TH1D *) fs_543->Get("MpPip");
    TH1D * h543_1 = (TH1D *) fs_543->Get("MpPim");
    TH1D * h543_2 = (TH1D *) fs_543->Get("MPiPi");
    TH1D * h543_4 = (TH1D *) fs_543->Get("Q2");

    TH2D * h543_5 = (TH2D *) fs_543->Get("W_vs_q2");
    TH2D * h543_5b = (TH2D *) fs_543->Get("W_vs_q2_photon");

    TH1D * h543_6a = (TH1D *) fs_543->Get("W_photon");
    TH1D * h543_6b = (TH1D *) fs_543->Get("W");
    
    TH1D * h543_7a = (TH1D *) fs_543->Get("W(03<q2<05)");
    TH1D * h543_7b = (TH1D *) fs_543->Get("W(03<q2<05)_photon");
    
    //error histogram 
    TH1D * herror_6a = new TH1D("W_photon_error", "W_photon_error", n_bin, 0, 5.5); 
    TH1D * herror_6b = new TH1D("W_error", "W_error", n_bin, 0, 5.5);

    float error_6a = 0;
    float error_6b = 0;
    for(int i=1; i<n_bin + 1; i++){
        error_6a = hdata_6a->GetBinContent(i);
        error_6b = hdata_7b->GetBinContent(i);
        
        error_6a = sqrt(error_6a);
        error_6b = sqrt(error_6b);
        
        herror_6a->SetBinContent(i, error_6a);
        herror_6b->SetBinContent(i, error_6b);
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
    
    h519_7a->Add(h523_7a,1);
    h519_7a->Add(h525_7a,1);
    h519_7a->Add(h543_7a,1);

    h519_7b->Add(h523_7b,1);
    h519_7b->Add(h525_7b,1);
    h519_7b->Add(h543_7b,1);
    //acceptance = simulation / original
    hsimu_0->Divide(h519_0);
    hsimu_1->Divide(h519_1);
    hsimu_2->Divide(h519_2);
    hsimu_4->Divide(h519_4);
    hsimu_5->Divide(h519_5);
    hsimu_5b->Divide(h519_5b);
    hsimu_6a->Divide(h519_6a);
    hsimu_6b->Divide(h519_6b);
    hsimu_7a->Divide(h519_7a);
    hsimu_7b->Divide(h519_7b);

    hsimu_0->SetDirectory(fs);
    hsimu_1->SetDirectory(fs);
    hsimu_2->SetDirectory(fs);
    hsimu_4->SetDirectory(fs);
    hsimu_5->SetDirectory(fs);
    hsimu_5b->SetDirectory(fs);
    hsimu_6a->SetDirectory(fs);
    hsimu_6b->SetDirectory(fs);
    hsimu_7a->SetDirectory(fs);
    hsimu_7b->SetDirectory(fs);

    hsimu_0->GetYaxis()->SetTitle("Acceptance");
    hsimu_1->GetYaxis()->SetTitle("Acceptance");
    hsimu_2->GetYaxis()->SetTitle("Acceptance");
    hsimu_4->GetYaxis()->SetTitle("Acceptance");
    hsimu_6a->GetYaxis()->SetTitle("Acceptance");
    hsimu_6b->GetYaxis()->SetTitle("Acceptance");
    hsimu_7a->GetYaxis()->SetTitle("Acceptance");
    hsimu_7b->GetYaxis()->SetTitle("Acceptance");

    hsimu_5->GetYaxis()->SetTitle("Q^2");
    hsimu_5b->GetYaxis()->SetTitle("Q^2");
    hsimu_4->GetXaxis()->SetTitle("W");
    hsimu_6a->GetXaxis()->SetTitle("W");
    hsimu_6b->GetXaxis()->SetTitle("W");
    
    hsimu_7a->GetXaxis()->SetTitle("W");
    hsimu_7b->GetXaxis()->SetTitle("W");
    
    hsimu_4->GetXaxis()->SetRangeUser(0.0, 1.0);

    hsimu_0->SetName("acceptance_invmass_Ppip");
    hsimu_1->SetName("acceptance_invmass_Ppim");
    hsimu_2->SetName("acceptance_invmass_pippim");
    hsimu_4->SetName("acceptance_Q2");
    hsimu_5->SetName("acceptance_WvsQ2");
    hsimu_5b->SetName("acceptance_WvsQ2_photon");
    hsimu_6a->SetName("acceptance_W_photon");
    hsimu_6b->SetName("acceptance_W");
    
    hsimu_7a->SetName("acceptance_W(02<q2<05)");
    hsimu_7b->SetName("acceptance_W(02<q2<05)_photon");
//cross section
    hdata_0->Divide(hsimu_0);
    hdata_1->Divide(hsimu_1);
    hdata_2->Divide(hsimu_2);
    hdata_4->Divide(hsimu_4);
    hdata_5->Divide(hsimu_5);
    hdata_5b->Divide(hsimu_5b);
    hdata_6a->Divide(hsimu_6a);
    hdata_6b->Divide(hsimu_6b);
    hdata_7a->Divide(hsimu_7a);
    hdata_7b->Divide(hsimu_7b);

    herror_6a->Divide(hsimu_6a);
    herror_6b->Divide(hsimu_6b);

    hdata_4->GetXaxis()->SetRangeUser(0.0, 1.0);

    hdata_0->Scale(1/convert);
    hdata_1->Scale(1/convert);
    hdata_2->Scale(1/convert);
    hdata_4->Scale(1/convert);
    hdata_5->Scale(1/convert);    
    hdata_5b->Scale(1/convert);
    hdata_6a->Scale(1/convert);
    hdata_6b->Scale(1/convert);
    hdata_7a->Scale(1/convert);
    hdata_7b->Scale(1/convert);

    herror_6a->Scale(1/convert);
    herror_6b->Scale(1/convert);
//binwidth
    hdata_0->Scale(n_bin/(3.5-1.0));
    hdata_1->Scale(n_bin/(3.5-1.0));
    hdata_2->Scale(n_bin/2.5);
    hdata_4->Scale(n_bin/4.5);
    hdata_5->Scale(n_bin/4.5*n_bin/5.5);
    hdata_5b->Scale(n_bin/4.5*n_bin/5.5);
    hdata_6a->Scale(n_bin/5.5);
    hdata_6b->Scale(n_bin/5.5);
    hdata_7a->Scale(n_bin/5.5);
    hdata_7b->Scale(n_bin/5.5);

    herror_6a->Scale(n_bin/5.5);
    herror_6b->Scale(n_bin/5.5);
//unit from cm2 to ub
    hdata_0->Scale(1e30);
    hdata_1->Scale(1e30);
    hdata_2->Scale(1e30);
    hdata_4->Scale(1e30); 
    hdata_5->Scale(1e30);
    hdata_5b->Scale(1e30);
    hdata_6a->Scale(1e30);
    hdata_6b->Scale(1e30);    
    hdata_7a->Scale(1e30);
    hdata_7b->Scale(1e30);

    herror_6a->Scale(1e30);
    herror_6b->Scale(1e30);

    for(int i=1; i<n_bin + 1; i++){
        error_6a = herror_6a->GetBinContent(i);
        error_6b = herror_6b->GetBinContent(i);
        
        hdata_6a->SetBinError(i, error_6a);
        hdata_6b->SetBinError(i, error_6b);
    }

    hdata_0->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_1->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_2->GetYaxis()->SetTitle("#frac{d#sigma}{dM} (#mub / GeV)");
    hdata_4->GetYaxis()->SetTitle("#frac{d#sigma}{dQ^2} (#mu b / GeV^2)");
    hdata_5->GetYaxis()->SetTitle("Q^2 (GeV^2)");
    hdata_5->GetXaxis()->SetTitle("W (GeV)");
    hdata_5b->GetYaxis()->SetTitle("Q^2 (GeV^2)");
    hdata_5b->GetXaxis()->SetTitle("W (GeV)");
    hdata_6a->GetYaxis()->SetTitle("#sigma_v (#mu b)");
    hdata_6b->GetYaxis()->SetTitle("#frac{d#sigma}{dW} (#mu b/ GeV)");
    hdata_7a->GetYaxis()->SetTitle("#frac{d#sigma}{dW} (#mu b/ GeV)");
    hdata_7b->GetYaxis()->SetTitle("#sigma_v (#mu b)");

    hdata_5->SetName("WvsQ2");
    hdata_5b->SetName("WvsQ2_photon");

    hdata_6a->SetName("cs_photon_W");
    hdata_6b->SetName("cs_elecprod_W");
    hdata_7a->SetName("cs_elecprod_W(0.3<q2<0.5)");
    hdata_7b->SetName("cs_photon_W(0.3<q2<0.5)");

//write
    hdata_0->SetDirectory(fs);
    hdata_1->SetDirectory(fs);
    hdata_2->SetDirectory(fs);
    hdata_4->SetDirectory(fs);
    hdata_5->SetDirectory(fs);   
    hdata_5b->SetDirectory(fs);
    hdata_6a->SetDirectory(fs);   
    hdata_6b->SetDirectory(fs);
    hdata_7a->SetDirectory(fs); 
    hdata_7b->SetDirectory(fs);

    fs->Write();
     
     return 1;
}
