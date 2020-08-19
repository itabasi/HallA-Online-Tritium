#include <iostream>
#include <TVector.h>
//#include <fstream>
//#include <string>
//#include <math.h>

using namespace std;
void mmass_simc_recon(){

  //  string rname = "/home/itabashi/SIMC/rootfiles/test.root";
  //  string rname = "/home/itabashi/SIMC/rootfiles/1H_Lam_0611.root";
  //    string rname = "./simc_root/3H_kaon_wb.root";
  string rname = "../rootfiles/simc/1H_kaon_new.root";

    //    bool nnL_mode =false;
    bool nnL_mode =true;
  
  TChain* T =new TChain("SNT");
  T->Add(Form("%s",rname.c_str()));

	  

  int RevID,LevID;
  float Rp_gen,Rth_gen,Rph_gen,Rx_gen,Ry_gen,Rz_gen;
  float Lp_gen,Lth_gen,Lph_gen,Lx_gen,Ly_gen,Lz_gen;
  
  float Rp_rec,Rth_rec,Rph_rec,Rx_rec,Ry_rec,Rz_rec;
  float Lp_rec,Lth_rec,Lph_rec,Lx_rec,Ly_rec,Lz_rec;    
  float Pm,nu;
  float uq_x,uq_y,uq_z,q;
    T->SetBranchAddress("Pm",&Pm);
    T->SetBranchAddress("uq_x_gen",&uq_x);
    T->SetBranchAddress("uq_y_gen",&uq_y);
    T->SetBranchAddress("uq_z_gen",&uq_z);
    T->SetBranchAddress("q_gen",&q);
    T->SetBranchAddress("nu",&nu);
    T->SetBranchAddress("Rp_gen",&Rp_gen);
    T->SetBranchAddress("Rth_gen",&Rth_gen);
    T->SetBranchAddress("Rph_gen",&Rth_gen);
    //    T->SetBranchAddress("h_xptar",&Rth_gen);
    //    T->SetBranchAddress("h_yptar",&Rph_gen);
    T->SetBranchAddress("Lp_gen",&Lp_gen);
    T->SetBranchAddress("Lth_gen",&Lth_gen);
    T->SetBranchAddress("Lph_gen",&Lth_gen);
    //    T->SetBranchAddress("e_xptar",&Lth_gen);
    //    T->SetBranchAddress("e_yptar",&Lph_gen);
    
    

    int ENum=T->GetEntries();


    double Rpx,Rpy,Rpz,Lpx,Lpy,Lpz;
    double Ek,Ee_,Ee,Bp;
    double Me   =  0.510998928e-3;
    double MT   =  2.808921112;
    double MnnL =  2.9948138266;
    double MK   =  0.493677;
    double Mp   =  0.938272046;
    double ML   =  1.115683;
    double mmass,mmass_2;
    string ofname="./simc_root/recon/Lam_Tmass.root";
    TFile* ofr =new TFile(ofname.c_str(),"recreate");
    TTree* Tnew=new TTree("T","T");
    Tnew= T->CloneTree(0);
    Tnew->Branch("mmass",&mmass,"mmass/D");
    Tnew->Branch("mmass_2",&mmass_2,"mmass/D");
    //    if(ENum>100000)ENum=100000;
    cout<<"ENum : "<<ENum<<endl;
    
    for(int i=0;i<ENum;i++){

      mmass=-10000.;
      mmass_2=-10000.;
      
      T->GetEntry(i);

      Ee  = 4.318;
      Lp_gen *= 1./1000.;
      Rp_gen *= 1./1000.;
      
      Ee_ = sqrt(Lp_gen * Lp_gen + Me*Me);
      Ek  = sqrt(Rp_gen * Rp_gen + MK*MK);
      Bp=sqrt(Ee*Ee-Me*Me);

      Rpz = Rp_gen / sqrt(1.0*1.0 + Rth_gen*Rth_gen + Rph_gen*Rph_gen );
      Rpx = Rpz * Rth_gen;
      Rpy = Rpz * Rph_gen;
      Lpz = Lp_gen / sqrt(1.0*1.0 + Lth_gen*Lth_gen + Lph_gen*Lph_gen );
      Lpx = Lpz * Lth_gen;
      Lpy = Lpz * Lph_gen;      

      TVector3 RP,LP,BP;
      BP.SetXYZ(0.0,0.0,Bp);
      RP.SetXYZ(Rpx,Rpy,Rpz);      
      LP.SetXYZ(Lpx,Lpy,Lpz);      

      //cout<<" Rpx "<<Rpx<<" Rpy "<<Rpy<<" Rpz "<<Rpz<<" Lpx "<<Lpx<<" Lpy "<<Lpy<<" Lpz "<<Lpz<<endl;
      
      RP.RotateX( -13.2/180.*3.14);
      LP.RotateX( +13.2/180.*3.14);


           
      if(nnL_mode){
      mmass = sqrt( (Ee + MT  - Ee_ - Ek )*(Ee + MT - Ee_ - Ek )
		    -(BP - LP - RP )*(BP - LP - RP)  );

      mmass = (mmass - MnnL)*1000.;
      //      mmass_2 = sqrt( ( nu +MT-Ek )*( nu + MT -Ek ) - Pm*Pm);
      mmass_2 = (nu + Mp - Ek)* (nu + Mp - Ek) ;
      mmass_2 = (mmass_2 -MnnL)*1000.;      
      
      }else{
      mmass = sqrt( (Ee + Mp  - Ee_ - Ek )*(Ee + Mp - Ee_ - Ek )
		    -(BP - LP - RP )*(BP - LP - RP)  );

      mmass = (mmass - ML)*1000.;

      mmass_2 = sqrt( ( nu +Mp-Ek )*( nu + Mp -Ek ) - Pm*Pm);
      mmass_2 = (mmass_2 -ML)*1000.;
      
      }
      //      cout<<"Ee "<<Ee<<" Ee_ "<<Ee_<<" Ek "<<Ek<<" MT "<<MT<<" a "<<sqrt((Ee + MT  - Ee_ - Ek )*(Ee + MT - Ee_ - Ek ))<<" b "<<sqrt((BP - LP - RP )*(BP - LP - RP) )<<" mmass "<<mmass<<endl; 
      
      //      cout<<"mass "<<mmass<<" Bp "<<Bp<<" Rp "<<Rp_gen<<" Lp "<<Lp_gen<<" Rth "<<Rth_gen*180./3.14<<" Lth "<<Lth_gen*180./3.14<<" Rph "<<Rph_gen*180./3.14<<" Lph "<<Lph_gen*180./3.14<<endl;
      Tnew->Fill();
    }

    Tnew->Write();

    cout<<"New RootFiles : "<<ofname<<endl;
  
}


