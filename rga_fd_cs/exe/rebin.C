#include "Lcore.h"


int main(){
    
    int ngroup =2; 
    
    TString filename1 = "output3_all_rebin.root";

    TFile * fs1 = new TFile(filename1.Data(), "RECREATE");

    //data
    TFile * fs_data = new TFile("/work/clas12/users/gjwei/tmp/2pion/build_new/output3_all.root", "r");
    TH1D * hdata_0 = (TH1D *) fs_data->Get("invmass_Ppip_twoPi");
    TH1D * hdata_1 = (TH1D *) fs_data->Get("invmass_Ppim_twoPi");
    TH1D * hdata_2 = (TH1D *) fs_data->Get("invmass_pippim_twoPi");
    TH1D * hdata_4 = (TH1D *) fs_data->Get("Q2_twoPi");
    TH2D * hdata_5 = (TH2D *) fs_data->Get("W_vs_q2_twoPi");
    TH1D * hdata_6a = (TH1D *) fs_data->Get("W_q2_15_twoPi");
    TH1D * hdata_6b = (TH1D *) fs_data->Get("W_q2_2_twoPi");
    TH1D * hdata_6c = (TH1D *) fs_data->Get("W_q2_23_twoPi");
    TH1D * hdata_6d = (TH1D *) fs_data->Get("W_q2_34_twoPi");
    TH1D * hdata_7a = (TH1D *) fs_data->Get("W_q2_15_photon_twoPi");
    TH1D * hdata_7b = (TH1D *) fs_data->Get("W_q2_2_photon_twoPi");
    TH1D * hdata_7c = (TH1D *) fs_data->Get("W_q2_23_photon_twoPi");
    TH1D * hdata_7d = (TH1D *) fs_data->Get("W_q2_34_photon_twoPi");

    TH1D * hdata_0_0 = (TH1D *) (hdata_0->Rebin(ngroup,"invmass_Ppip_twoPi"));
    TH1D * hdata_1_1 = (TH1D *) (hdata_1->Rebin(ngroup,"invmass_Ppim_twoPi"));
    TH1D * hdata_2_2 = (TH1D *) (hdata_2->Rebin(ngroup,"invmass_pippim_twoPi"));
    TH1D * hdata_4_4 = (TH1D *) (hdata_4->Rebin(ngroup,"Q2_twoPi"));
    TH2D * hdata_5_5 = (TH2D *) (hdata_5->Rebin2D(ngroup,ngroup,"W_vs_q2_twoPi"));
    TH1D * hdata_6a_6a = (TH1D *) (hdata_6a->Rebin(ngroup,"W_q2_15_twoPi"));
    TH1D * hdata_6b_6b = (TH1D *) (hdata_6b->Rebin(ngroup,"W_q2_2_twoPi"));
    TH1D * hdata_6c_6c = (TH1D *) (hdata_6c->Rebin(ngroup,"W_q2_23_twoPi"));
    TH1D * hdata_6d_6d = (TH1D *) (hdata_6d->Rebin(ngroup,"W_q2_34_twoPi"));
    TH1D * hdata_7a_7a = (TH1D *) (hdata_7a->Rebin(ngroup,"W_q2_15_photon_twoPi"));
    TH1D * hdata_7b_7b = (TH1D *) (hdata_7b->Rebin(ngroup,"W_q2_2_photon_twoPi"));
    TH1D * hdata_7c_7c = (TH1D *) (hdata_7c->Rebin(ngroup,"W_q2_23_photon_twoPi"));
    TH1D * hdata_7d_7d = (TH1D *) (hdata_7d->Rebin(ngroup,"W_q2_34_photon_twoPi"));

    hdata_0_0->SetDirectory(fs1);
    hdata_1_1->SetDirectory(fs1);
    hdata_2_2->SetDirectory(fs1);
    hdata_4_4->SetDirectory(fs1);
    hdata_5_5->SetDirectory(fs1);
    hdata_6a_6a->SetDirectory(fs1);
    hdata_6b_6b->SetDirectory(fs1);
    hdata_6c_6c->SetDirectory(fs1);
    hdata_6d_6d->SetDirectory(fs1);
    hdata_7a_7a->SetDirectory(fs1);
    hdata_7b_7b->SetDirectory(fs1);
    hdata_7c_7c->SetDirectory(fs1);
    hdata_7d_7d->SetDirectory(fs1);
    fs1->Write();
     
    TString filename2 = "output_new_simu_rebin.root";

    TFile * fs2 = new TFile(filename2.Data(), "RECREATE");
    //simu
    TFile * fs_simu = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_script/before_rebin/output_new_simu.root", "r");
    
    TH1D * hsimu_0 = (TH1D *) fs_simu->Get("invmass_Ppip_twoPi");
    TH1D * hsimu_1 = (TH1D *) fs_simu->Get("invmass_Ppim_twoPi");
    TH1D * hsimu_2 = (TH1D *) fs_simu->Get("invmass_pippim_twoPi");
    TH1D * hsimu_4 = (TH1D *) fs_simu->Get("Q2_twoPi");
    TH2D * hsimu_5 = (TH2D *) fs_simu->Get("W_vs_q2_twoPi");
    TH1D * hsimu_6a = (TH1D *) fs_simu->Get("W_q2_15_twoPi");
    TH1D * hsimu_6b = (TH1D *) fs_simu->Get("W_q2_2_twoPi");
    TH1D * hsimu_6c = (TH1D *) fs_simu->Get("W_q2_23_twoPi");
    TH1D * hsimu_6d = (TH1D *) fs_simu->Get("W_q2_34_twoPi");
    TH1D * hsimu_7a = (TH1D *) fs_simu->Get("W_q2_15_photon_twoPi");
    TH1D * hsimu_7b = (TH1D *) fs_simu->Get("W_q2_2_photon_twoPi");
    TH1D * hsimu_7c = (TH1D *) fs_simu->Get("W_q2_23_photon_twoPi");
    TH1D * hsimu_7d = (TH1D *) fs_simu->Get("W_q2_34_photon_twoPi"); 

    TH1D * hsimu_0_0 = (TH1D *) (hsimu_0->Rebin(ngroup,"invmass_Ppip_twoPi"));
    TH1D * hsimu_1_1 = (TH1D *) (hsimu_1->Rebin(ngroup,"invmass_Ppim_twoPi"));
    TH1D * hsimu_2_2 = (TH1D *) (hsimu_2->Rebin(ngroup,"invmass_pippim_twoPi"));
    TH1D * hsimu_4_4 = (TH1D *) (hsimu_4->Rebin(ngroup,"Q2_twoPi"));
    TH2D * hsimu_5_5 = (TH2D *) (hsimu_5->Rebin2D(ngroup,ngroup,"W_vs_q2_twoPi"));
    TH1D * hsimu_6a_6a = (TH1D *) (hsimu_6a->Rebin(ngroup,"W_q2_15_twoPi"));
    TH1D * hsimu_6b_6b = (TH1D *) (hsimu_6b->Rebin(ngroup,"W_q2_2_twoPi"));
    TH1D * hsimu_6c_6c = (TH1D *) (hsimu_6c->Rebin(ngroup,"W_q2_23_twoPi"));
    TH1D * hsimu_6d_6d = (TH1D *) (hsimu_6d->Rebin(ngroup,"W_q2_34_twoPi"));
    TH1D * hsimu_7a_7a = (TH1D *) (hsimu_7a->Rebin(ngroup,"W_q2_15_photon_twoPi"));
    TH1D * hsimu_7b_7b = (TH1D *) (hsimu_7b->Rebin(ngroup,"W_q2_2_photon_twoPi"));
    TH1D * hsimu_7c_7c = (TH1D *) (hsimu_7c->Rebin(ngroup,"W_q2_23_photon_twoPi"));
    TH1D * hsimu_7d_7d = (TH1D *) (hsimu_7d->Rebin(ngroup,"W_q2_34_photon_twoPi"));

    hsimu_0_0->SetDirectory(fs2);
    hsimu_1_1->SetDirectory(fs2);
    hsimu_2_2->SetDirectory(fs2);
    hsimu_4_4->SetDirectory(fs2);
    hsimu_5_5->SetDirectory(fs2);
    hsimu_6a_6a->SetDirectory(fs2);
    hsimu_6b_6b->SetDirectory(fs2);
    hsimu_6c_6c->SetDirectory(fs2);
    hsimu_6d_6d->SetDirectory(fs2);
    hsimu_7a_7a->SetDirectory(fs2);
    hsimu_7b_7b->SetDirectory(fs2);
    hsimu_7c_7c->SetDirectory(fs2);
    hsimu_7d_7d->SetDirectory(fs2);
    fs2->Write();

    TString filename3 = "original_574_rebin.root";

    TFile * fs3 = new TFile(filename3.Data(), "RECREATE");
    
    TFile * fs_519 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_script/new_lund/original_574.root", "r");

    TH1D * h519_0 = (TH1D *) fs_519->Get("MpPip");
    TH1D * h519_1 = (TH1D *) fs_519->Get("MpPim");
    TH1D * h519_2 = (TH1D *) fs_519->Get("MPiPi");
    TH1D * h519_4 = (TH1D *) fs_519->Get("Q2");
    TH2D * h519_5 = (TH2D *) fs_519->Get("W_vs_q2");
    TH1D * h519_6a = (TH1D *) fs_519->Get("W(q2<1.5)");
    TH1D * h519_6b = (TH1D *) fs_519->Get("W(q2<2)");
    TH1D * h519_6c = (TH1D *) fs_519->Get("W(2<q2<3)");
    TH1D * h519_6d = (TH1D *) fs_519->Get("W(3<q2<4)");
    TH2D * h519_7a = (TH2D *) fs_519->Get("W(q2<1.5)_photon");
    TH1D * h519_7b = (TH1D *) fs_519->Get("W(q2<2)_photon");
    TH1D * h519_7c = (TH1D *) fs_519->Get("W(2<q2<3)_photon");
    TH1D * h519_7d = (TH1D *) fs_519->Get("W(3<q2<4)_photon");
    
    TH1D * h519_0_0 = (TH1D *) (h519_0->Rebin(ngroup,"MpPip"));
    TH1D * h519_1_1 = (TH1D *) (h519_1->Rebin(ngroup,"MpPim"));
    TH1D * h519_2_2 = (TH1D *) (h519_2->Rebin(ngroup,"MPiPi"));
    TH1D * h519_4_4 = (TH1D *) (h519_4->Rebin(ngroup,"Q2"));
    TH2D * h519_5_5 = (TH2D *) (h519_5->Rebin2D(ngroup,ngroup,"W_vs_q2"));
    TH1D * h519_6a_6a = (TH1D *) (h519_6a->Rebin(ngroup,"W(q2<1.5)"));
    TH1D * h519_6b_6b = (TH1D *) (h519_6b->Rebin(ngroup,"W(q2<2)"));
    TH1D * h519_6c_6c = (TH1D *) (h519_6c->Rebin(ngroup,"W(2<q2<3)"));
    TH1D * h519_6d_6d = (TH1D *) (h519_6d->Rebin(ngroup,"W(3<q2<4)"));
    TH1D * h519_7a_7a = (TH1D *) (h519_7a->Rebin(ngroup,"W(q2<1.5)_photon"));
    TH1D * h519_7b_7b = (TH1D *) (h519_7b->Rebin(ngroup,"W(q2<2)_photon"));
    TH1D * h519_7c_7c = (TH1D *) (h519_7c->Rebin(ngroup,"W(2<q2<3)_photon"));
    TH1D * h519_7d_7d = (TH1D *) (h519_7d->Rebin(ngroup,"W(3<q2<4)_photon"));

    h519_0_0->SetDirectory(fs3);
    h519_1_1->SetDirectory(fs3);
    h519_2_2->SetDirectory(fs3);
    h519_4_4->SetDirectory(fs3);
    h519_5_5->SetDirectory(fs3);
    h519_6a_6a->SetDirectory(fs3);
    h519_6b_6b->SetDirectory(fs3);
    h519_6c_6c->SetDirectory(fs3);
    h519_6d_6d->SetDirectory(fs3);
    h519_7a_7a->SetDirectory(fs3);
    h519_7b_7b->SetDirectory(fs3);
    h519_7c_7c->SetDirectory(fs3);
    h519_7d_7d->SetDirectory(fs3);
    fs3->Write();

    TString filename4 = "original_575_rebin.root";

    TFile * fs4 = new TFile(filename4.Data(), "RECREATE");
    
    TFile * fs_523 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_script/new_lund/original_575.root", "r");
    TH1D * h523_0 = (TH1D *) fs_523->Get("MpPip");
    TH1D * h523_1 = (TH1D *) fs_523->Get("MpPim");
    TH1D * h523_2 = (TH1D *) fs_523->Get("MPiPi");
    TH1D * h523_4 = (TH1D *) fs_523->Get("Q2");
    TH2D * h523_5 = (TH2D *) fs_523->Get("W_vs_q2");
    TH1D * h523_6a = (TH1D *) fs_523->Get("W(q2<1.5)");
    TH1D * h523_6b = (TH1D *) fs_523->Get("W(q2<2)");
    TH1D * h523_6c = (TH1D *) fs_523->Get("W(2<q2<3)");
    TH1D * h523_6d = (TH1D *) fs_523->Get("W(3<q2<4)");
    TH2D * h523_7a = (TH2D *) fs_523->Get("W(q2<1.5)_photon");
    TH1D * h523_7b = (TH1D *) fs_523->Get("W(q2<2)_photon");
    TH1D * h523_7c = (TH1D *) fs_523->Get("W(2<q2<3)_photon");
    TH1D * h523_7d = (TH1D *) fs_523->Get("W(3<q2<4)_photon");
    
    TH1D * h523_0_0 = (TH1D *) (h523_0->Rebin(ngroup,"MpPip"));
    TH1D * h523_1_1 = (TH1D *) (h523_1->Rebin(ngroup,"MpPim"));
    TH1D * h523_2_2 = (TH1D *) (h523_2->Rebin(ngroup,"MPiPi"));
    TH1D * h523_4_4 = (TH1D *) (h523_4->Rebin(ngroup,"Q2"));
    TH2D * h523_5_5 = (TH2D *) (h523_5->Rebin2D(ngroup,ngroup,"W_vs_q2"));
    TH1D * h523_6a_6a = (TH1D *) (h523_6a->Rebin(ngroup,"W(q2<1.5)"));
    TH1D * h523_6b_6b = (TH1D *) (h523_6b->Rebin(ngroup,"W(q2<2)"));
    TH1D * h523_6c_6c = (TH1D *) (h523_6c->Rebin(ngroup,"W(2<q2<3)"));
    TH1D * h523_6d_6d = (TH1D *) (h523_6d->Rebin(ngroup,"W(3<q2<4)"));
    TH1D * h523_7a_7a = (TH1D *) (h523_7a->Rebin(ngroup,"W(q2<1.5)_photon"));
    TH1D * h523_7b_7b = (TH1D *) (h523_7b->Rebin(ngroup,"W(q2<2)_photon"));
    TH1D * h523_7c_7c = (TH1D *) (h523_7c->Rebin(ngroup,"W(2<q2<3)_photon"));
    TH1D * h523_7d_7d = (TH1D *) (h523_7d->Rebin(ngroup,"W(3<q2<4)_photon"));

    h523_0_0->SetDirectory(fs4);
    h523_1_1->SetDirectory(fs4);
    h523_2_2->SetDirectory(fs4);
    h523_4_4->SetDirectory(fs4);
    h523_5_5->SetDirectory(fs4);
    h523_6a_6a->SetDirectory(fs4);
    h523_6b_6b->SetDirectory(fs4);
    h523_6c_6c->SetDirectory(fs4);
    h523_6d_6d->SetDirectory(fs4);
    h523_7a_7a->SetDirectory(fs4);
    h523_7b_7b->SetDirectory(fs4);
    h523_7c_7c->SetDirectory(fs4);
    h523_7d_7d->SetDirectory(fs4);
    fs4->Write();

    TString filename5 = "original_576_rebin.root";

    TFile * fs5 = new TFile(filename5.Data(), "RECREATE");
    
    TFile * fs_525 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_script/new_lund/original_576.root", "r");
    TH1D * h525_0 = (TH1D *) fs_525->Get("MpPip");
    TH1D * h525_1 = (TH1D *) fs_525->Get("MpPim");
    TH1D * h525_2 = (TH1D *) fs_525->Get("MPiPi");
    TH1D * h525_4 = (TH1D *) fs_525->Get("Q2");
    TH2D * h525_5 = (TH2D *) fs_525->Get("W_vs_q2");
    TH1D * h525_6a = (TH1D *) fs_525->Get("W(q2<1.5)");
    TH1D * h525_6b = (TH1D *) fs_525->Get("W(q2<2)");
    TH1D * h525_6c = (TH1D *) fs_525->Get("W(2<q2<3)");
    TH1D * h525_6d = (TH1D *) fs_525->Get("W(3<q2<4)");
    TH2D * h525_7a = (TH2D *) fs_525->Get("W(q2<1.5)_photon");
    TH1D * h525_7b = (TH1D *) fs_525->Get("W(q2<2)_photon");
    TH1D * h525_7c = (TH1D *) fs_525->Get("W(2<q2<3)_photon");
    TH1D * h525_7d = (TH1D *) fs_525->Get("W(3<q2<4)_photon");
    
    TH1D * h525_0_0 = (TH1D *) (h525_0->Rebin(ngroup,"MpPip"));
    TH1D * h525_1_1 = (TH1D *) (h525_1->Rebin(ngroup,"MpPim"));
    TH1D * h525_2_2 = (TH1D *) (h525_2->Rebin(ngroup,"MPiPi"));
    TH1D * h525_4_4 = (TH1D *) (h525_4->Rebin(ngroup,"Q2"));
    TH2D * h525_5_5 = (TH2D *) (h525_5->Rebin2D(ngroup,ngroup,"W_vs_q2"));
    TH1D * h525_6a_6a = (TH1D *) (h525_6a->Rebin(ngroup,"W(q2<1.5)"));
    TH1D * h525_6b_6b = (TH1D *) (h525_6b->Rebin(ngroup,"W(q2<2)"));
    TH1D * h525_6c_6c = (TH1D *) (h525_6c->Rebin(ngroup,"W(2<q2<3)"));
    TH1D * h525_6d_6d = (TH1D *) (h525_6d->Rebin(ngroup,"W(3<q2<4)"));
    TH1D * h525_7a_7a = (TH1D *) (h525_7a->Rebin(ngroup,"W(q2<1.5)_photon"));
    TH1D * h525_7b_7b = (TH1D *) (h525_7b->Rebin(ngroup,"W(q2<2)_photon"));
    TH1D * h525_7c_7c = (TH1D *) (h525_7c->Rebin(ngroup,"W(2<q2<3)_photon"));
    TH1D * h525_7d_7d = (TH1D *) (h525_7d->Rebin(ngroup,"W(3<q2<4)_photon"));

    h525_0_0->SetDirectory(fs5);
    h525_1_1->SetDirectory(fs5);
    h525_2_2->SetDirectory(fs5);
    h525_4_4->SetDirectory(fs5);
    h525_5_5->SetDirectory(fs5);
    h525_6a_6a->SetDirectory(fs5);
    h525_6b_6b->SetDirectory(fs5);
    h525_6c_6c->SetDirectory(fs5);
    h525_6d_6d->SetDirectory(fs5);
    h525_7a_7a->SetDirectory(fs5);
    h525_7b_7b->SetDirectory(fs5);
    h525_7c_7c->SetDirectory(fs5);
    h525_7d_7d->SetDirectory(fs5);
    fs5->Write();

    TString filename6 = "original_577_rebin.root";

    TFile * fs6 = new TFile(filename6.Data(), "RECREATE");
    
    TFile * fs_543 = new TFile("/work/clas12/users/gjwei/acceptance/code_original/rga_script/new_lund/original_577.root", "r");
    TH1D * h543_0 = (TH1D *) fs_543->Get("MpPip");
    TH1D * h543_1 = (TH1D *) fs_543->Get("MpPim");
    TH1D * h543_2 = (TH1D *) fs_543->Get("MPiPi");
    TH1D * h543_4 = (TH1D *) fs_543->Get("Q2");
    TH2D * h543_5 = (TH2D *) fs_543->Get("W_vs_q2");
    TH1D * h543_6a = (TH1D *) fs_543->Get("W(q2<1.5)");
    TH1D * h543_6b = (TH1D *) fs_543->Get("W(q2<2)");
    TH1D * h543_6c = (TH1D *) fs_543->Get("W(2<q2<3)");
    TH1D * h543_6d = (TH1D *) fs_543->Get("W(3<q2<4)");
    TH2D * h543_7a = (TH2D *) fs_543->Get("W(q2<1.5)_photon");
    TH1D * h543_7b = (TH1D *) fs_543->Get("W(q2<2)_photon");
    TH1D * h543_7c = (TH1D *) fs_543->Get("W(2<q2<3)_photon");
    TH1D * h543_7d = (TH1D *) fs_543->Get("W(3<q2<4)_photon");
    
    TH1D * h543_0_0 = (TH1D *) (h543_0->Rebin(ngroup,"MpPip"));
    TH1D * h543_1_1 = (TH1D *) (h543_1->Rebin(ngroup,"MpPim"));
    TH1D * h543_2_2 = (TH1D *) (h543_2->Rebin(ngroup,"MPiPi"));
    TH1D * h543_4_4 = (TH1D *) (h543_4->Rebin(ngroup,"Q2"));
    TH2D * h543_5_5 = (TH2D *) (h543_5->Rebin2D(ngroup,ngroup,"W_vs_q2"));
    TH1D * h543_6a_6a = (TH1D *) (h543_6a->Rebin(ngroup,"W(q2<1.5)"));
    TH1D * h543_6b_6b = (TH1D *) (h543_6b->Rebin(ngroup,"W(q2<2)"));
    TH1D * h543_6c_6c = (TH1D *) (h543_6c->Rebin(ngroup,"W(2<q2<3)"));
    TH1D * h543_6d_6d = (TH1D *) (h543_6d->Rebin(ngroup,"W(3<q2<4)"));
    TH1D * h543_7a_7a = (TH1D *) (h543_7a->Rebin(ngroup,"W(q2<1.5)_photon"));
    TH1D * h543_7b_7b = (TH1D *) (h543_7b->Rebin(ngroup,"W(q2<2)_photon"));
    TH1D * h543_7c_7c = (TH1D *) (h543_7c->Rebin(ngroup,"W(2<q2<3)_photon"));
    TH1D * h543_7d_7d = (TH1D *) (h543_7d->Rebin(ngroup,"W(3<q2<4)_photon"));

    h543_0_0->SetDirectory(fs6);
    h543_1_1->SetDirectory(fs6);
    h543_2_2->SetDirectory(fs6);
    h543_4_4->SetDirectory(fs6);
    h543_5_5->SetDirectory(fs6);
    h543_6a_6a->SetDirectory(fs6);
    h543_6b_6b->SetDirectory(fs6);
    h543_6c_6c->SetDirectory(fs6);
    h543_6d_6d->SetDirectory(fs6);
    h543_7a_7a->SetDirectory(fs6);
    h543_7b_7b->SetDirectory(fs6);
    h543_7c_7c->SetDirectory(fs6);
    h543_7d_7d->SetDirectory(fs6);
    fs6->Write();

   return 1;
}
