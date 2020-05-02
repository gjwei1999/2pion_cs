#include "Lcore.h"

int main(){

    gStyle->SetPalette(53);

    //data
    TFile * fs_data = new TFile("cs_rga_eft.root", "r"); 
      
    TH1D * hdata = (TH1D *) fs_data->Get("Q2_twoPi");

    TFile * fs_gen = new TFile("cs_generator.root", "r");

    TH1D * hgen = (TH1D *) fs_gen->Get("Q2");

    hdata->SetLineColor(2); hgen->SetLineColor(4);

    TLegend * leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg4->AddEntry(hdata, "From data", "l");
    leg4->AddEntry(hgen, "From generator", "l");
    TCanvas *c4 = new TCanvas("c4","",800,600);
    hdata->GetXaxis()->SetRangeUser(0.0,1.0);
    hdata->GetYaxis()->SetRangeUser(0.0,0.4);
    hdata->Draw();
    hgen->Draw("same");
    leg4->Draw("same");


return 1;
}




  
