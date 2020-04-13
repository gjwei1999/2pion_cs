#include "Lcore.h"

int main(){
    gStyle->SetPalette(53);

    TString filename = "rga_cs_generator.root";

    TFile * fs = new TFile(filename.Data(), "RECREATE");
    
    //original from generator
    TFile * fs_519 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_492.root", "r");
    TH1D * h519_0 = (TH1D *) fs_519->Get("MpPip");
    TH1D * h519_1 = (TH1D *) fs_519->Get("MpPim");
    TH1D * h519_2 = (TH1D *) fs_519->Get("MPiPi");
    TH1D * h519_3 = (TH1D *) fs_519->Get("MissingSQ");
    TH1D * h519_4 = (TH1D *) fs_519->Get("Q2");

    TFile * fs_523 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_493.root", "r");
    TH1D * h523_0 = (TH1D *) fs_523->Get("MpPip");
    TH1D * h523_1 = (TH1D *) fs_523->Get("MpPim");
    TH1D * h523_2 = (TH1D *) fs_523->Get("MPiPi");
    TH1D * h523_3 = (TH1D *) fs_523->Get("MissingSQ");
    TH1D * h523_4 = (TH1D *) fs_523->Get("Q2");

    TFile * fs_525 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_498.root", "r");
    TH1D * h525_0 = (TH1D *) fs_525->Get("MpPip");
    TH1D * h525_1 = (TH1D *) fs_525->Get("MpPim");
    TH1D * h525_2 = (TH1D *) fs_525->Get("MPiPi");
    TH1D * h525_3 = (TH1D *) fs_525->Get("MissingSQ");
    TH1D * h525_4 = (TH1D *) fs_525->Get("Q2");

    TFile * fs_543 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_ori_lund/original_499.root", "r");
    TH1D * h543_0 = (TH1D *) fs_543->Get("MpPip");
    TH1D * h543_1 = (TH1D *) fs_543->Get("MpPim");
    TH1D * h543_2 = (TH1D *) fs_543->Get("MPiPi");
    TH1D * h543_3 = (TH1D *) fs_543->Get("MissingSQ");
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

    h519_3->Add(h523_3,1);
    h519_3->Add(h525_3,1);
    h519_3->Add(h543_3,1);

    h519_4->Add(h523_4,1);
    h519_4->Add(h525_4,1);
    h519_4->Add(h543_4,1);
    
    float f0 = 2e7*(3.5-1.0)/500;
    float f1 = 2e7*(3.5-1.0)/500;
    float f2 = 2e7*(2.5)/500;
    float f4 = 2e7/192;
    
    h519_0->Scale(1/f0);
    h519_1->Scale(1/f1);
    h519_2->Scale(1/f2);
    h519_4->Scale(1/f4);

    h519_0->SetDirectory(fs);
    h519_1->SetDirectory(fs);
    h519_2->SetDirectory(fs);
    h519_4->SetDirectory(fs);
    
    h519_0->GetXaxis()->SetTitle("M(p#pi^{+}) (GeV)");
    h519_0->GetYaxis()->SetTitle("#frac{d#sigma}{dM}(#mu b /GeV)");
    
    h519_1->GetXaxis()->SetTitle("M(p#pi^{-}) (GeV)");
    h519_1->GetYaxis()->SetTitle("#frac{d#sigma}{dM}(#mu b /GeV)");
    
    h519_2->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) (GeV)");
    h519_2->GetYaxis()->SetTitle("#frac{d#sigma}{dM}(#mu b /GeV)");
    
    h519_4->GetXaxis()->SetTitle("Q^2 (GeV)");
    h519_4->GetYaxis()->SetTitle("#frac{#sigma}{dQ^2} (#mu b/GeV^2)");
    
    fs->Write();
    
    return 1;

}
