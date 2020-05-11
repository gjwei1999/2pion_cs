/****************************************/
/*                                      */
/*       Created by Jiawei Guo          */
/*           Duke University            */
/*                                      */
/****************************************/

#include "Lcore.h"

int main(const int argc, const char * argv[]){
    
    if (argc < 2){
    std::cout << "./projection <rebin ngroup> " << endl;
    return 0;
  }

    int ngroup = atoi(argv[1]);
    float n_bin=500 / ngroup;

    TString filename = "projection.root";

    TFile * fs = new TFile(filename.Data(), "RECREATE");

    TFile * fs_data = new TFile("/work/clas12/users/gjwei/rga_fd_cs/build/cs_rga1.root", "r");
    
    TH2D * hdata_5a = (TH2D *) fs_data->Get("WvsQ2");
    TH2D * hdata_5b = (TH2D *) fs_data->Get("WvsQ2_photon");

    int i1, i2;
    i1 = hdata_5a->GetYaxis()->FindBin(2.0);
    i2 = hdata_5b->GetYaxis()->FindBin(3.0);
    TH1D * hprojection_5a = (TH1D *) hdata_5a->ProjectionX("W(2<q2<3)", i1, i2, "e");
    

    i1 = hdata_5a->GetYaxis()->FindBin(2.0);
    i2 = hdata_5b->GetYaxis()->FindBin(3.0);
    TH1D * hprojection_5b = (TH1D *) hdata_5b->ProjectionX("W(2<q2<3)_photon", i1, i2, "e");

    hprojection_5a->Scale(12.0/n_bin);
    hprojection_5b->Scale(12.0/n_bin);

    hprojection_5a->SetDirectory(fs);
    hprojection_5b->SetDirectory(fs);

    hprojection_5b->GetYaxis()->SetRangeUser(0.0, 3.0);   
    fs->Write();    
return 1;
}
