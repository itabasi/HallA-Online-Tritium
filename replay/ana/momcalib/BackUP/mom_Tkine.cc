


void mom_Tkine(){
  string rname="../rootfiles/momcalib/test.root";

  TChain * T = new TChain("T");
  T->Add(rname.c_str());
  

  
  int ENum = T->GetEntries();
  cout<<"ENum : "<<ENum<<endl;


  



}
