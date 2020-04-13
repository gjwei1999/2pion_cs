#include "Lcore.h"

int main(){

    gStyle->SetPalette(53);

//    TString filename = "crosssection.root";

  //  TFile * fs = new TFile(filename.Data(), "RECREATE");

    //data
    TFile * fs_data = new TFile("/work/clas12/users/gjwei/multi_thr_clas12_ana/build_new/output3_new.root", "r");

    TH1D * hdata_0 = (TH1D *) fs_data->Get("invmass_Ppip_twoPi");
    TH1D * hdata_1 = (TH1D *) fs_data->Get("invmass_Ppim_twoPi");
    TH1D * hdata_2 = (TH1D *) fs_data->Get("invmass_pippim_twoPi");
    TH1D * hdata_4 = (TH1D *) fs_data->Get("Q2_twoPi");

    //simulation    
    TFile * fs_simu = new TFile("/work/clas12/users/gjwei/acceptance/rga/output_rga_simu.root", "r");
    TH1D * hsimu_0 = (TH1D *) fs_simu->Get("invmass_Ppip_twoPi");                     
    TH1D * hsimu_1 = (TH1D *) fs_simu->Get("invmass_Ppim_twoPi");                            
    TH1D * hsimu_2 = (TH1D *) fs_simu->Get("invmass_pippim_twoPi");                              
    TH1D * hsimu_4 = (TH1D *) fs_simu->Get("Q2_twoPi");

//original from generator
    TFile * fs_519 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_492.root", "r");     
    TH1D * h519_0 = (TH1D *) fs_519->Get("MpPip");
    TH1D * h519_1 = (TH1D *) fs_519->Get("MpPim");     
    TH1D * h519_2 = (TH1D *) fs_519->Get("MPiPi");                 
    TH1D * h519_4 = (TH1D *) fs_519->Get("Q2");        

    TFile * fs_523 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_493.root", "r");
    TH1D * h523_0 = (TH1D *) fs_523->Get("MpPip");
    TH1D * h523_1 = (TH1D *) fs_523->Get("MpPim");                                  
    TH1D * h523_2 = (TH1D *) fs_523->Get("MPiPi");
    TH1D * h523_4 = (TH1D *) fs_523->Get("Q2");                                 
    
    TFile * fs_525 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_498.root", "r");     TH1D * h525_0 = (TH1D *) fs_525->Get("MpPip");
    TH1D * h525_1 = (TH1D *) fs_525->Get("MpPim");
    TH1D * h525_2 = (TH1D *) fs_525->Get("MPiPi");
    TH1D * h525_4 = (TH1D *) fs_525->Get("Q2");

    TFile * fs_543 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_499.root", "r");
    TH1D * h543_0 = (TH1D *) fs_543->Get("MpPip");
    TH1D * h543_1 = (TH1D *) fs_543->Get("MpPim");
    TH1D * h543_2 = (TH1D *) fs_543->Get("MPiPi");
    TH1D * h543_4 = (TH1D *) fs_543->Get("Q2");
                                                                
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
//
    h519_0->SetLineColor(2); hsimu_0->SetLineColor(4);
    h519_1->SetLineColor(2); hsimu_1->SetLineColor(4);
    h519_2->SetLineColor(2); hsimu_2->SetLineColor(4);
    h519_4->SetLineColor(2); hsimu_4->SetLineColor(4);

    TLegend * leg0 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg0->AddEntry(h519_0, "Generated", "l");
    leg0->AddEntry(hsimu_0, "Simulation", "l");    
    TCanvas *c0 = new TCanvas("c0","",800,600);
    h519_0->Draw();
    hsimu_0->Draw("same");
    leg0->Draw("same");
    
    TLegend * leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(h519_1, "Generated", "l");
    leg1->AddEntry(hsimu_1, "Simulation", "l");
    TCanvas *c1 = new TCanvas("c1","",800,600);
    h519_1->Draw();
    hsimu_1->Draw("same");
    leg1->Draw("same");

    TLegend * leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->AddEntry(h519_2, "Generated", "l");
    leg2->AddEntry(hsimu_2, "Simulation", "l");
    TCanvas *c2 = new TCanvas("c2","",800,600);
    h519_2->Draw();
    hsimu_2->Draw("same");
    leg2->Draw("same");

    TLegend * leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg4->AddEntry(h519_4, "Generated", "l");
    leg4->AddEntry(hsimu_4, "Simulation", "l");
    TCanvas *c4 = new TCanvas("c4","",800,600);
    h519_4->Draw();
    hsimu_4->Draw("same");
    leg4->Draw("same");    
return 1;
}
