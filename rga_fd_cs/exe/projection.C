/****************************************/
/*                                      */
/*       Created by Jiawei Guo          */
/*           Duke University            */
/*                                      */
/****************************************/

#include "Lcore.h"

int main(){
    
    TString filename = "projection.root";

    TFile * fs = new TFile(filename.Data(), "RECREATE");

    TFile * fs_data = new TFile("/work/clas12/users/gjwei/rga_fd_cs/build/cs_rga.root", "r");
    
    TH2D * hdata_5 = (TH2D *) fs_data->Get("WvsQ2");
    
    int i1, i2;
    i1 = hdata_5->GetYaxis()->FindBin(2.0);
    i2 = hdata_5->GetYaxis()->FindBin(3.0);
    TH1D * hprojection = (TH1D *) hdata_5->ProjectionX("W(2<q2<3)", i1, i2, "e");
    
    hprojection->SetDirectory(fs);


    fs->Write();    
return 1;
}
