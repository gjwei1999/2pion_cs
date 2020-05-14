/****************************************/
/*                                      */
/*       Created by Jiawei Guo          */
/*           Duke University            */
/*                                      */
/****************************************/

#include "Lcore.h"

int main(){

    gStyle->SetPalette(53);

    //data
    TFile * fs_data = new TFile("cs_rga.root", "r"); 
      
    TH1D * hdata_0 = (TH1D *) fs_data->Get("cs_photon_W(2<q2<3)");
    TH1D * hdata_1 = (TH1D *) fs_data->Get("cs_photon_W(3<q2<4)");
    //TFile * fs_gen = new TFile("cs_generator.root", "r");

    //TH1D * hgen = (TH1D *) fs_gen->Get("Q2");

    hdata_0->SetLineColor(2); hdata_1->SetLineColor(4);

    TLegend * leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg4->AddEntry(hdata_0, "W(2<q2<3)", "l");
    leg4->AddEntry(hdata_1, "W(3<q2<4)", "l");
    TCanvas *c4 = new TCanvas("c4","",800,600);
    //hdata->GetXaxis()->SetRangeUser(0.0,1.0);
    //hdata->GetYaxis()->SetRangeUser(0.0,0.4);
    hdata_0->Draw();
    hdata_1->Draw("same");
    leg4->Draw("same");


return 1;
}




  
