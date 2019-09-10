
#include "Setting.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  int ch;
  extern char *optarg;
  string ifname("input.root");
  string ofname("output.root");

  while((ch=getopt(argc,argv,"hf:w:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input root :"<<ifname<<endl;
      break;
    case 'w':
      ofname = optarg;
      cout<<"output root :"<<ofname<<endl;
      break;
    case 'h':
      std::cout<<"-f (inputfile): input ROOT file name"<<std::endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  
  Setting *set = new Setting();
  set->Initialize();

 TApplication theApp("App", &argc, argv);

 //======== TChain ==================//

  TChain *oldtree = new TChain("T");
  //  oldtree->Add(Form("%s",ifname.c_str()));
  
 ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    oldtree->Add(runname.c_str());
    //   cout<<buf<<endl;
  }

 


  double L_tr_n, L_tr_x[20], L_tr_th[20], L_tr_p[20];
  double R_tr_n, R_tr_x[20], R_tr_th[20], R_tr_p[20];

  oldtree->SetBranchStatus("*",0);
  

 oldtree->SetBranchStatus("RTDC.F1FirstHit",1);
 oldtree->SetBranchStatus("R.s2.t_pads",1);
 oldtree->SetBranchStatus("R.s2.trpad",1);
 oldtree->SetBranchStatus("R.a1.asum_c",1);
 oldtree->SetBranchStatus("R.a2.asum_c",1);
  // path length//
 oldtree->SetBranchStatus("R.s2.trpath",1); 
 oldtree->SetBranchStatus("R.tr.pathl",1);  
 // Target positon information //
 oldtree->SetBranchStatus("R.tr.p",1);
 oldtree->SetBranchStatus("R.tr.vz",1);    

 //------ Left Arm ---------------//
 oldtree->SetBranchStatus("LTDC.F1FirstHit",1);
 oldtree->SetBranchStatus("L.s2.t_pads",1);
 oldtree->SetBranchStatus("L.s2.trpad",1);
  // path length//
 oldtree->SetBranchStatus("L.s2.trpath",1); 
 oldtree->SetBranchStatus("L.tr.pathl",1);   
 oldtree->SetBranchStatus("L.tr.p",1);
 oldtree->SetBranchStatus("L.tr.vz",1);    

  

  TFile *ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  TTree *newtree = oldtree->CloneTree(1);
  int ENum = oldtree->GetEntries(0);

  
  for(int n=0;n<ENum;n++){
    oldtree->GetEntry(n);

    int NLtr = (int)L_tr_n;
    int NRtr = (int)R_tr_n;
    bool LOK = false;
    bool ROK = false;

    for(int lt=0;lt<NLtr;lt++){
      if(   L_tr_th[lt]<0.17*L_tr_x[lt]+0.050
         && L_tr_th[lt]>0.17*L_tr_x[lt]-0.060
         && L_tr_p[lt]>1. && L_tr_p[lt]<10.) LOK = true;
    }

    for(int rt=0;rt<NRtr;rt++){
      if(   R_tr_th[rt]<0.17*R_tr_x[rt]+0.050
         && R_tr_th[rt]>0.17*R_tr_x[rt]-0.060
         && R_tr_p[rt]>1. && R_tr_p[rt]<10. ) ROK = true;
    }

    if( LOK && ROK ){
      //      newtree->Fill();
    }
    newtree->Fill(); //optics run
    if(n % 100000 == 0){ cout<<n<<" / "<<ENum<<endl; }

  } // for ENum

  
  newtree->AutoSave();
  delete ofp;

  gSystem->Exit(1);
  theApp.Run();
  return 0;

} 



