

void check_ole_root(){

  string old_root="/data3/root/tritium_111132.root";
  string ole_root="/data4/root/tritium_111132.root";

  TFile * ifr_old =new TFile(old_root.c_str());
  TFile * ifr_ole =new TFile(ole_root.c_str());

  TChain* Told=new TChain("T");
  TChain* Tole=new TChain("T");

  double Rz_ole,Rz_old;
  
  Told->Add(old_root.c_str());
  Tole->Add(ole_root.c_str());

  Told->SetBranchStatus("*",0);
  Told->SetBranchStatus("R.tr.vz",1);  
  Told->SetBranchAddress("R.tr.vz",&Rz_old);

  Tole->SetBranchStatus("*",0);
  Tole->SetBranchStatus("R.tr.vz",1);  
  Tole->SetBranchAddress("R.tr.vz",&Rz_ole);
  //  Tole->SetBranchStatus("fEvtHdr.fEvtNum",1);  
  Tole->SetBranchAddress("R.tr.vz",&Rz_ole);
  int ENum_old=Told->GetEntries();
  int ENum_ole=Tole->GetEntries();

  
  string ofname ="./check.root";
  TFile* ofs=new TFile(ofname.c_str(),"recreate");  
  TTree* tnew = new TTree("T","check ROOT");
  tnew->Branch("Rz_old",&Rz_old,"Rz_old/D");
  tnew->Branch("Rz_ole",&Rz_ole,"Rz_ole/D");
  

  cout<<"ENum_old "<<ENum_old<<endl;
  for(int i=0;i<ENum_old;i++){
    Rz_old=-1000.;
    Told->GetEntry(i);
    tnew->Fill();
    cout<<"i "<<i<<endl;
  }
  cout<<"ENum_ole "<<ENum_ole<<endl;
  for(int i=0;i<ENum_ole;i++){
    Rz_ole=-1000.;
    Tole->GetEntry(i);
    tnew->Fill();
  }


  tnew->Write();
}



