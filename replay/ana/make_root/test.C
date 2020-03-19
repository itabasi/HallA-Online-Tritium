//#include "THaEvent.C"
void test(){
  //gSystem->Load("/home/gaio/analyzer-1.6.3/libHallA.so.1.6");
  //gSystem->Load("/home/gaio/analyzer-1.6.3/libHallA.so.1.6");
  TFile *file = new TFile("/data3/root_ole/root2/tritium_111200_1.root","readonly");
  TTree *T=(TTree*)file->Get("T");
  unsigned int evnum;
  T->SetBranchAddress("fEvtHdr.fEvtNum", &evnum);

  int ENum = T->GetEntries();
  ENum=10;
  for(int n=0;n<ENum;n++){
    T->GetEntry(n);
    cout<<n<<" "<<evnum<<endl;
  }

}
