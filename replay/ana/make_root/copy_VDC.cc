
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

  TChain *oldtree = new TChain("T");
  oldtree->Add(Form("%s",ifname.c_str()));


  oldtree->SetBranchStatus("*",0);
  
  oldtree->SetBranchStatus("fEvtHdr.fRun"   ,1);
  oldtree->SetBranchStatus("fEvtHdr.fEvtNum" ,1);
  oldtree->SetBranchStatus("HALLA_p"        ,1);
  oldtree->SetBranchStatus("DR.evtypebits"        ,1);
  //================//
  //====== VDC =====//
  //===============//

  //========== RHRS VDC ==========//
  
  oldtree->SetBranchStatus("R.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.nhit"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);    
  oldtree->SetBranchStatus("R.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.dist"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.ddist"          ,1);    
  oldtree->SetBranchStatus("R.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.trdist"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.nclust"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.nclust"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.nclust"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.nclust"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.clsiz"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.clsiz"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.clsiz"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.clsiz"          ,1);   
  oldtree->SetBranchStatus("R.vdc.u1.clpivot"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.clpivot"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.clpivot"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.clpivot"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.clpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.clpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.clpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.clpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.slope"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.slope"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.slope"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.slope"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.lslope"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.lslope"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.lslope"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.lslope"          ,1);
  /*
  oldtree->SetBranchStatus("R.vdc.u1.t0"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.t0"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.t0"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.t0"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.sigsl"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.sigsl"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.sigsl"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.sigsl"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.sigpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.sigpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.sigpos"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.sigpos"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.sigt0"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.sigt0"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.sigt0"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.sigt0"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.clchi2"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.clchi2"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.clchi2"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.clchi2"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.clndof"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.clndof"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.clndof"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.clndof"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.cltcor"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.cltcor"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.cltcor"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.cltcor"          ,1);  
  */
  
  //========== LHRS VDC ===============//

  oldtree->SetBranchStatus("L.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.nhit"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.rawtime"          ,1);    
  oldtree->SetBranchStatus("L.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.time"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.dist"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.ddist"          ,1);    
  oldtree->SetBranchStatus("L.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.trdist"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.nclust"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.nclust"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.nclust"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.nclust"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.nclust"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.clsiz"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.clsiz"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.clsiz"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.clsiz"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.clsiz"          ,1);   
  oldtree->SetBranchStatus("L.vdc.u1.clpivot"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.clpivot"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.clpivot"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.clpivot"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.clpivot"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.clpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.clpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.clpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.clpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.clpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.u1.slope"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.slope"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.slope"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.slope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.slope"          ,1);
  oldtree->SetBranchStatus("L.vdc.u1.lslope"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.lslope"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.lslope"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.lslope"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.lslope"          ,1);
  /*
  oldtree->SetBranchStatus("L.vdc.u1.t0"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.t0"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.t0"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.t0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.t0"          ,1);
  oldtree->SetBranchStatus("L.vdc.u1.sigsl"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.sigsl"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.sigsl"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.sigsl"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.sigsl"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.sigpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.sigpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.sigpos"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.sigpos"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.sigpos"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.sigt0"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.sigt0"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.sigt0"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.sigt0"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.sigt0"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.clchi2"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.clchi2"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.clchi2"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.clchi2"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.clchi2"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.clndof"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.clndof"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.clndof"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.clndof"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.clndof"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.cltcor"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.cltcor"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.cltcor"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.cltcor"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.cltcor"          ,1);  
  */
  
  
  TFile *ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  TTree *newtree = oldtree->CloneTree();
  newtree->Write();
  int ENum = oldtree->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

  
  //  newtree->AutoSave();
  ofp->Close();
    delete ofp;

  gSystem->Exit(1);
  theApp.Run();
  return 0;

} 



