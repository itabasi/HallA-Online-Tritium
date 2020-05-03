#include "Setting.h"
#include "TObject.h"
#include "TBranchElement.h"
using namespace std;
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
extern double s2f1_off(int i,string ARM, int KINE);

class EvtHdr : public TObject
{

public :
    ULong64_t EvtTime;
    UInt_t EvtNum;
    Int_t EvtType;
    Int_t EvtLen;
    Int_t Helicity;
    Int_t TargetPol;
    Int_t Run;

  //  EvtHdr(){};
  //  virtual ~EvtHdr(){};
  //  ClassDef(EvtHdr,1);
};


//ClassInp(EvtHdr)

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

  TTree* oldtree2 =new TTree(ifname.c_str(),"T");

  double L_tr_n, L_tr_x[20], L_tr_th[20], L_tr_p[20];
  double R_tr_n, R_tr_x[20], R_tr_th[20], R_tr_p[20];
  double R_s2_trpad[20],L_s2_trpad[20];
  double R_tr_pathl[20],L_tr_pathl[20],LF1[200],RF1[200],R_s2_trpath[20],L_s2_trpath[20];
  double ct[10],rtof[10],ltof[10];
  int nev,runnum;
  UInt_t Nev;
  int Nrun;


  

  /*
  typedef struct {
    ULong64_t EvtTime;
    UInt_t EvtNum;
    Int_t EvtType;
    Int_t EvtLen;
    Int_t Helicity;
    Int_t TargetPol;
    Int_t Run;
  } fEvtHdr;
  fEvtHdr test;
  */

  //  oldtree->SetBranchStatus("*"        ,0);  
  //  oldtree->SetBranchStatus("fEvtHdr.fEvtNum"        ,1);



    //  EvtHdr* test = new EvtHdr();
  
  oldtree->SetBranchStatus("*",0);
  oldtree->SetBranchStatus("fEvtHdr.fEvtNum"        ,1);  
  oldtree->SetBranchStatus("HALLA_p"        ,1);
  oldtree->SetBranchStatus("HALLA_dpp"        ,1);

  //========== Beam Raster =====================//
  //  oldtree->SetBranchStatus("Rrb.Raster2.target.x" ,1);
  //  oldtree->SetBranchStatus("Rrb.Raster2.target.y" ,1);
  oldtree->SetBranchStatus("Rrb.Raster2.rawcur.x" ,1);
  oldtree->SetBranchStatus("Rrb.Raster2.rawcur.y" ,1);
  //  oldtree->SetBranchStatus("FbusRrb.Raster2.rawcur.x" ,1);
  //  oldtree->SetBranchStatus("FbusRrb.Raster2.rawcur.y" ,1);
  //  oldtree->SetBranchStatus("Lrb.Raster2.target.x" ,1);
  //  oldtree->SetBranchStatus("Lrb.Raster2.target.y" ,1);
  oldtree->SetBranchStatus("Lrb.Raster2.rawcur.x" ,1);
  oldtree->SetBranchStatus("Lrb.Raster2.rawcur.y" ,1);
  //  oldtree->SetBranchStatus("FbusLrb.Raster2.rawcur.x" ,1);
  //  oldtree->SetBranchStatus("FbusLrb.Raster2.rawcur.y" ,1);
  oldtree->SetBranchStatus("rbax"   ,1);
  oldtree->SetBranchStatus("rbay"   ,1);
  oldtree->SetBranchStatus("rbbx"   ,1);
  oldtree->SetBranchStatus("rbby"   ,1);
  oldtree->SetBranchStatus("rbx"    ,1);
  oldtree->SetBranchStatus("rby"    ,1);
  oldtree->SetBranchStatus("bpmaws" ,1);
  oldtree->SetBranchStatus("bpmbws" ,1);
  oldtree->SetBranchStatus("DL.evtype"       ,1);
  oldtree->SetBranchStatus("L.s0.la"         ,1);
  oldtree->SetBranchStatus("L.s0.la_c"       ,1);
  oldtree->SetBranchStatus("L.s0.la_p"       ,1);
  oldtree->SetBranchStatus("L.s0.lt"         ,1);
  oldtree->SetBranchStatus("L.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra"         ,1);
  oldtree->SetBranchStatus("L.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s0.rt"         ,1);
  oldtree->SetBranchStatus("L.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("L.s0.time"       ,1);
  oldtree->SetBranchStatus("L.s0.dedx"       ,1);
  oldtree->SetBranchStatus("L.s0.trdy"       ,1);
  oldtree->SetBranchStatus("L.s0.troff"      ,1);
  oldtree->SetBranchStatus("L.s0.trpad"      ,1); 
  oldtree->SetBranchStatus("L.s0.trpath"     ,1);
  oldtree->SetBranchStatus("L.s0.trx"        ,1);
  oldtree->SetBranchStatus("L.s0.try"        ,1);
  oldtree->SetBranchStatus("L.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s2.la"         ,1);
  oldtree->SetBranchStatus("L.s2.la_c"       ,1);
  oldtree->SetBranchStatus("L.s2.la_p"       ,1);
  oldtree->SetBranchStatus("L.s2.lt"         ,1);
  oldtree->SetBranchStatus("L.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra"         ,1);
  oldtree->SetBranchStatus("L.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s2.rt"         ,1);
  oldtree->SetBranchStatus("L.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.t_pads"     ,1); 
  oldtree->SetBranchStatus("L.s2.time"       ,1);
  oldtree->SetBranchStatus("L.s2.dedx"       ,1);
  oldtree->SetBranchStatus("L.s2.trdx"       ,1);
  //  oldtree->SetBranchStatus("L.s2.troff"      ,1);
  oldtree->SetBranchStatus("L.s2.trpad"      ,1); oldtree->SetBranchAddress("L.s2.trpad",L_s2_trpad);
  oldtree->SetBranchStatus("L.s2.trpath"     ,1); oldtree->SetBranchAddress("L.s2.trpath",L_s2_trpath);
  oldtree->SetBranchStatus("L.s2.trx"        ,1);
  oldtree->SetBranchStatus("L.s2.try"        ,1);
  oldtree->SetBranchStatus("L.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("L.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.u1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.u2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.v1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.v2.nhit"          ,1);  
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
  oldtree->SetBranchStatus("L.cer.a"         ,1);
  oldtree->SetBranchStatus("L.cer.a_c"       ,1);
  oldtree->SetBranchStatus("L.cer.a_p"       ,1);
  oldtree->SetBranchStatus("L.cer.t"         ,1);
  oldtree->SetBranchStatus("L.cer.t_c"       ,1);
  oldtree->SetBranchStatus("L.cer.asum_p"    ,1);
  oldtree->SetBranchStatus("L.cer.asum_c"    ,1);
  oldtree->SetBranchStatus("L.cer.trpath"    ,1);
  oldtree->SetBranchStatus("L.cer.trx"       ,1);
  oldtree->SetBranchStatus("L.cer.try"       ,1);
  oldtree->SetBranchStatus("L.cer.nhits"     ,1);
  oldtree->SetBranchStatus("L.cer.peak"      ,1);
  oldtree->SetBranchStatus("L.cer.t_fadc"    ,1);
  oldtree->SetBranchStatus("L.cer.tc_fadc"   ,1);
  //  oldtree->SetBranchStatus("L.ps.asum_p"     ,1);
  //  oldtree->SetBranchStatus("L.ps.asum_c"     ,1);
  //  oldtree->SetBranchStatus("L.sh.asum_p"     ,1);
  //  oldtree->SetBranchStatus("L.sh.asum_c"     ,1);
  oldtree->SetBranchStatus("L.tr.n"          ,1); oldtree->SetBranchAddress("L.tr.n" ,&L_tr_n );
  oldtree->SetBranchStatus("L.tr.flag"       ,1);
  oldtree->SetBranchStatus("L.tr.ndof"       ,1);
  oldtree->SetBranchStatus("L.tr.chi2"       ,1);
  oldtree->SetBranchStatus("L.tr.beta"       ,1);
  oldtree->SetBranchStatus("L.tr.d_x"        ,1);
  oldtree->SetBranchStatus("L.tr.d_y"        ,1);
  oldtree->SetBranchStatus("L.tr.d_th"       ,1);
  oldtree->SetBranchStatus("L.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.r_x"        ,1);
  oldtree->SetBranchStatus("L.tr.r_y"        ,1);
  oldtree->SetBranchStatus("L.tr.r_th"       ,1);
  oldtree->SetBranchStatus("L.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.x"          ,1);  oldtree->SetBranchAddress("L.tr.x" ,L_tr_x );
  oldtree->SetBranchStatus("L.tr.y"          ,1);
  oldtree->SetBranchStatus("L.tr.th"         ,1);  oldtree->SetBranchAddress("L.tr.th",L_tr_th);
  oldtree->SetBranchStatus("L.tr.ph"         ,1);
  oldtree->SetBranchStatus("L.tr.time"       ,1);
  oldtree->SetBranchStatus("L.tr.p"          ,1);  oldtree->SetBranchAddress("L.tr.p" ,L_tr_p );
  oldtree->SetBranchStatus("L.tr.pathl"      ,1);  oldtree->SetBranchAddress("L.tr.pathl" ,L_tr_pathl);
  oldtree->SetBranchStatus("L.tr.px"         ,1);
  oldtree->SetBranchStatus("L.tr.py"         ,1);
  oldtree->SetBranchStatus("L.tr.pz"         ,1);
  oldtree->SetBranchStatus("L.tr.tg_dp"      ,1);
  oldtree->SetBranchStatus("L.tr.tg_y"       ,1);
  oldtree->SetBranchStatus("L.tr.tg_th"      ,1);
  oldtree->SetBranchStatus("L.tr.tg_ph"      ,1);
  oldtree->SetBranchStatus("L.tr.vx"         ,1);
  oldtree->SetBranchStatus("L.tr.vy"         ,1);
  oldtree->SetBranchStatus("L.tr.vz"         ,1);
  oldtree->SetBranchStatus("DL.ltRFtime"     ,1);
  oldtree->SetBranchStatus("LTDC.F1FirstHit" ,1); oldtree->SetBranchAddress("LTDC.F1FirstHit" ,LF1);
  oldtree->SetBranchStatus("Ndata.DR.T1"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T2"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T3"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T4"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T5"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T6"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T7"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T8"     ,1);
  oldtree->SetBranchStatus("DR.T1"     ,1);
  oldtree->SetBranchStatus("DR.T2"     ,1);
  oldtree->SetBranchStatus("DR.T3"     ,1);
  oldtree->SetBranchStatus("DR.T4"     ,1);
  oldtree->SetBranchStatus("DR.T5"     ,1);
  oldtree->SetBranchStatus("DR.T6"     ,1);
  oldtree->SetBranchStatus("DR.T7"     ,1);
  oldtree->SetBranchStatus("DR.T8"     ,1);
  oldtree->SetBranchStatus("DR.evtype"       ,1);
  oldtree->SetBranchStatus("DR.evtypebits"   ,1);  
  oldtree->SetBranchStatus("R.s0.la"         ,1);
  oldtree->SetBranchStatus("R.s0.la_c"       ,1);
  oldtree->SetBranchStatus("R.s0.la_p"       ,1);
  oldtree->SetBranchStatus("R.s0.lt"         ,1);
  oldtree->SetBranchStatus("R.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra"         ,1);
  oldtree->SetBranchStatus("R.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s0.rt"         ,1);
  oldtree->SetBranchStatus("R.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s0.time"       ,1);
  oldtree->SetBranchStatus("R.s0.dedx"       ,1);
  oldtree->SetBranchStatus("R.s0.trdy"       ,1);
  oldtree->SetBranchStatus("R.s0.troff"      ,1);
  oldtree->SetBranchStatus("R.s0.trpad"      ,1);
  oldtree->SetBranchStatus("R.s0.trpath"     ,1);
  oldtree->SetBranchStatus("R.s0.trx"        ,1);
  oldtree->SetBranchStatus("R.s0.try"        ,1);
  oldtree->SetBranchStatus("R.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s2.la"         ,1);
  oldtree->SetBranchStatus("R.s2.la_c"       ,1);
  oldtree->SetBranchStatus("R.s2.la_p"       ,1);
  oldtree->SetBranchStatus("R.s2.lt"         ,1);
  oldtree->SetBranchStatus("R.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra"         ,1);
  oldtree->SetBranchStatus("R.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s2.rt"         ,1);
  oldtree->SetBranchStatus("R.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s2.time"       ,1);
  oldtree->SetBranchStatus("R.s2.dedx"       ,1);
  oldtree->SetBranchStatus("R.s2.trdx"       ,1);
  //  oldtree->SetBranchStatus("R.s2.troff"      ,1);
  oldtree->SetBranchStatus("R.s2.trpad"      ,1);  oldtree->SetBranchAddress("R.s2.trpad",R_s2_trpad);
  oldtree->SetBranchStatus("R.s2.trpath"     ,1);  oldtree->SetBranchAddress("R.s2.trpath",R_s2_trpath);
  oldtree->SetBranchStatus("R.s2.trx"        ,1);
  oldtree->SetBranchStatus("R.s2.try"        ,1);
  oldtree->SetBranchStatus("R.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("R.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.u1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.u2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.v1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.v2.nhit"          ,1);  
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
  oldtree->SetBranchStatus("R.a1.a"          ,1);
  oldtree->SetBranchStatus("R.a1.a_c"        ,1);
  oldtree->SetBranchStatus("R.a1.a_p"        ,1);
  oldtree->SetBranchStatus("R.a1.t"          ,1);
  oldtree->SetBranchStatus("R.a1.t_c"        ,1);
  oldtree->SetBranchStatus("R.a1.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a1.asum_c"     ,1);
  oldtree->SetBranchStatus("R.a1.trpath"     ,1);
  oldtree->SetBranchStatus("R.a1.trx"        ,1);
  oldtree->SetBranchStatus("R.a1.try"        ,1);
  oldtree->SetBranchStatus("R.a1.nhits"      ,1);
  oldtree->SetBranchStatus("R.a1.peak"       ,1);
  oldtree->SetBranchStatus("R.a1.t_fadc"     ,1);
  oldtree->SetBranchStatus("R.a1.tc_fadc"    ,1);
  oldtree->SetBranchStatus("R.a1.nahit"      ,1);     
  oldtree->SetBranchStatus("R.a2.a"          ,1);
  oldtree->SetBranchStatus("R.a2.a_c"        ,1);
  oldtree->SetBranchStatus("R.a2.a_p"        ,1);
  oldtree->SetBranchStatus("R.a2.t"          ,1);
  oldtree->SetBranchStatus("R.a2.t_c"        ,1);
  oldtree->SetBranchStatus("R.a2.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a2.asum_c"     ,1);
  oldtree->SetBranchStatus("R.a2.trpath"     ,1);
  oldtree->SetBranchStatus("R.a2.trx"        ,1);
  oldtree->SetBranchStatus("R.a2.try"        ,1);
  oldtree->SetBranchStatus("R.a2.nhits"      ,1);
  oldtree->SetBranchStatus("R.a2.peak"       ,1);
  oldtree->SetBranchStatus("R.a2.t_fadc"     ,1);
  oldtree->SetBranchStatus("R.a2.tc_fadc"    ,1);
  oldtree->SetBranchStatus("R.a2.nahit"      ,1); 

  /*
  oldtree->SetBranchStatus("R.cer.a"         ,1);
  oldtree->SetBranchStatus("R.cer.a_c"       ,1);
  oldtree->SetBranchStatus("R.cer.a_p"       ,1);
  oldtree->SetBranchStatus("R.cer.t"         ,1);
  oldtree->SetBranchStatus("R.cer.t_c"       ,1);
  oldtree->SetBranchStatus("R.cer.asum_p"    ,1);
  oldtree->SetBranchStatus("R.cer.asum_c"    ,1);
  oldtree->SetBranchStatus("R.cer.trpath"    ,1);
  oldtree->SetBranchStatus("R.cer.trx"       ,1);
  oldtree->SetBranchStatus("R.cer.try"       ,1);
  oldtree->SetBranchStatus("R.cer.nhits"     ,1);
  oldtree->SetBranchStatus("R.cer.peak"      ,1);
  oldtree->SetBranchStatus("R.cer.t_fadc"    ,1);
  oldtree->SetBranchStatus("R.cer.tc_fadc"   ,1);
  oldtree->SetBranchStatus("R.ps.asum_p"     ,1);
  oldtree->SetBranchStatus("R.ps.asum_c"     ,1);
  oldtree->SetBranchStatus("R.sh.asum_p"     ,1);
  oldtree->SetBranchStatus("R.sh.asum_c"     ,1);  
  */
  oldtree->SetBranchStatus("R.tr.n"          ,1);  oldtree->SetBranchAddress("R.tr.n" ,&R_tr_n );
  oldtree->SetBranchStatus("R.tr.flag"       ,1);
  oldtree->SetBranchStatus("R.tr.ndof"       ,1);
  oldtree->SetBranchStatus("R.tr.chi2"       ,1);
  oldtree->SetBranchStatus("R.tr.beta"       ,1);
  oldtree->SetBranchStatus("R.tr.d_x"        ,1);
  oldtree->SetBranchStatus("R.tr.d_y"        ,1);
  oldtree->SetBranchStatus("R.tr.d_th"       ,1);
  oldtree->SetBranchStatus("R.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.r_x"        ,1);  oldtree->SetBranchAddress("R.tr.x" ,R_tr_x );
  oldtree->SetBranchStatus("R.tr.r_y"        ,1);
  oldtree->SetBranchStatus("R.tr.r_th"       ,1);  oldtree->SetBranchAddress("R.tr.th",R_tr_th);
  oldtree->SetBranchStatus("R.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.x"          ,1);
  oldtree->SetBranchStatus("R.tr.y"          ,1);  
  oldtree->SetBranchStatus("R.tr.th"         ,1);
  oldtree->SetBranchStatus("R.tr.ph"         ,1);
  oldtree->SetBranchStatus("R.tr.time"       ,1);
  oldtree->SetBranchStatus("R.tr.p"          ,1); oldtree->SetBranchAddress("R.tr.p" ,R_tr_p );
  oldtree->SetBranchStatus("R.tr.pathl"      ,1); oldtree->SetBranchAddress("R.tr.pathl" ,R_tr_pathl);
  oldtree->SetBranchStatus("R.tr.px"         ,1);
  oldtree->SetBranchStatus("R.tr.py"         ,1);
  oldtree->SetBranchStatus("R.tr.pz"         ,1);
  oldtree->SetBranchStatus("R.tr.tg_dp"      ,1);
  oldtree->SetBranchStatus("R.tr.tg_y"       ,1);
  oldtree->SetBranchStatus("R.tr.tg_th"      ,1);
  oldtree->SetBranchStatus("R.tr.tg_ph"      ,1);
  oldtree->SetBranchStatus("R.tr.vx"         ,1);
  oldtree->SetBranchStatus("R.tr.vy"         ,1);
  oldtree->SetBranchStatus("R.tr.vz"         ,1);
  oldtree->SetBranchStatus("DR.rtRFtime"     ,1);
  oldtree->SetBranchStatus("RTDC.F1FirstHit" ,1); oldtree->SetBranchAddress("RTDC.F1FirstHit" ,RF1);


  TBranchElement* evID = (TBranchElement*)oldtree->GetBranch("fEvtHdr.fEvtNum");


  //==== Get Run number =====//
  int nrun;
  string nrun_c,st;
  
    for(int i=0;i<(int)ifname.length();i++){
      st=ifname.substr(i,1);
      if(st=="1"){
	  nrun_c=ifname.substr(i,6);
	  nrun=atoi(nrun_c.c_str());
	  if(nrun>111111){
	    runnum=nrun;
	  }
      }
    }
    double tdc_time,coin_offset;
    int tdc_mode;
    const  double pathl_off=-470.5;
    if(runnum<=111368){
      tdc_time=0.056;// [ns]
      tdc_mode=1;
      coin_offset=pathl_off + 464.73; // H1 mode
    }else{
      tdc_time=0.058;// [ns]
      tdc_mode=2;
      coin_offset=pathl_off + 470.63; // H2 mode
    }


  


  TFile *ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  TTree *newtree = oldtree->CloneTree(0);


  oldtree->SetBranchStatus("fEvtHdr"   ,1);
  //  oldtree->SetBranchAddress("fEvtHdr" , &test);

  //  TBranch* evbranch= oldtree->GetBranch("fEvtHdr");
  //  oldtree->SetBranchAddress("fEvtHdr" , &test);
  //  oldtree->SetBranchAddress("fEvtHdr.fEvtNum" , &test->EvtNum);
    //    oldtree->SetBranchAddress("fEvtHdr" ,test.EvtTime:test.EvtNum:test.EvtType:test.EvtLen:test.Helicity:test.TargetPol:test.Run);



  newtree->Branch("ct",ct,"ct[10]/D");
  newtree->Branch("rtof",rtof,"rtof[10]/D");
  newtree->Branch("ltof",ltof,"ltof[10]/D");
  newtree->Branch("runnum",&runnum,"runnum/I");
  newtree->Branch("nev",&nev,"nev/I");
  //  newtree->Branch("fEvtHdr",&test,"EvtTime/D:EvtNum/I:EvtType/I:EvtLen/I:Helicity/I:TargetPol/I:Run/I");
  //  newtree->Branch("fEvtHdr",&test,"test.EvtTime/D:test.EvtNum/I:test.EvtType/I:test.EvtLen/I:test.Helicity/I:test.TargetPol/I:test.Run/I");



  /*

  typedef struct {
    ULong64_t EvtTime;
    UInt_t EvtNum;
    Int_t EvtType;
    Int_t EvtLen;
    Int_t Helicity;
    Int_t TargetPol;
    Int_t Run;
  } fEvtHdr;
  fEvtHdr test;

   */

  int ENum = oldtree->GetEntries();
  cout<<"ENum "<<ENum<<endl;

  for(int n=0;n<ENum;n++){
    for(int j=0;j<10;j++)ct[j]=-100.0;

    Nev=0;
    Nrun=0;

    /*
    test.EvtNum=0;
    test.EvtTime=0;
    test.EvtType=0;
    test.EvtLen=0;
    test.TargetPol=0;
    test.Run=0;
    */


    //    test->EvtNum=0;
    //    test->EvtTime=0;
    //    test->EvtType=0;
    //    test->EvtLen=0;
    //    test->TargetPol=0;
    //    test->Run=0;


    oldtree->GetEntry(n);
    evID->GetEntry(n);
    nev = evID->GetValue(0,0);
    //    nev=n;
    /*
    cout<<"n "<<n<<" nev "<<test.EvtNum<<" time "<<test.EvtTime<<" type "<<test.EvtType<<" len "<<test.EvtLen<<" helicity "<<test.Helicity<<" targetPol "<<test.TargetPol<<" run "<<test.Run<<endl;
    */

    //    cout<<"n "<<n<<" nev "<<test->EvtNum<<" time "<<test->EvtTime<<" type "<<test->EvtType<<" len "<<test->EvtLen<<" helicity "<<test->Helicity<<" targetPol "<<test->TargetPol<<" run "<<test->Run<<endl;


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
    int nct=0;
    for(int lt=0;lt<NLtr;lt++){
      for(int rt=0;rt<NRtr;rt++){
	int    Ls2pads=(int)L_s2_trpad[lt];
	int    Rs2pads=(int)R_s2_trpad[rt];
	double Rs2_off= s2f1_off(Rs2pads,"R",tdc_mode);
	double Ls2_off= s2f1_off(Ls2pads,"L",tdc_mode);
	/*
	double tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
	double tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
	double R_betaK=R_tr_p[rt]/sqrt(MK*MK + R_tr_p[rt]*R_tr_p[rt]);
	double L_tgt = tof_l + (L_tr_pathl[lt] + L_s2_trpath[lt])/c;
	double R_tgt = tof_r + (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaK/c;
	*/

	//======= F1 TDC Shift correction ===============//
	
	double RS2T = -RF1[16+Rs2pads]+RF1[9];
	double RS2B = -RF1[48+Rs2pads]+RF1[46];
	double LS2T = -LF1[Ls2pads]+LF1[30];
	double LS2B = -LF1[Ls2pads+48]+LF1[37];

	/*
	if(RS2T>pow(2,15))RS2T=RS2T -pow(2.0,16) + 350.;
	else if(RS2T<-pow(2,15))RS2T=RS2T + pow(2.0,16) -350.;
	if(RS2B>pow(2,15))RS2B=RS2B -pow(2.0,16) + 350.;
	else if(RS2B<-pow(2,15))RS2B=RS2B + pow(2.0,16) -350.;      
	if(LS2T>pow(2,15))LS2T=LS2T -pow(2.0,16) + 350.;
	else if(LS2T<-pow(2,15))LS2T=LS2T + pow(2.0,16) -350.;
	if(LS2B>pow(2,15))LS2B=LS2B -pow(2.0,16) + 350.;
	else if(LS2B<-pow(2,15))LS2B=LS2B + pow(2.0,16) -350.;
	*/
      
	double tof_r=( ( RS2T + RS2B + Rs2_off )/2.0 )*tdc_time;
	double tof_l=( ( LS2T + LS2B + Ls2_off )/2.0 )*tdc_time;
	double R_betaK=R_tr_p[rt]/sqrt(MK*MK + R_tr_p[rt]*R_tr_p[rt]);
	double L_tgt = tof_l + (L_tr_pathl[lt] + L_s2_trpath[lt])/c;
	double R_tgt = tof_r + (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaK/c;

	//====================================================//
	
	rtof[rt]=R_tgt-coin_offset;
	ltof[lt]=L_tgt;
	ct[nct] = - L_tgt + R_tgt- coin_offset;
	nct++;
      }
    }
    
    if( LOK && ROK ){
      newtree->Fill();
    }
    //newtree->Fill(); //optics run
    if(n % 100000 == 0){ cout<<n<<" / "<<ENum<<endl; }

  } // for ENum

  newtree->AutoSave();
  delete ofp;

  gSystem->Exit(1);
  theApp.Run();
  return 0;

}




// ####################################################
double s2f1_off(int i,string ARM, int KINE){
// ####################################################

  double RS2_offset[16],LS2_offset[16];
  if(KINE==2){

    double RS2_R_off[16]={-8361.42, -8395.25, -8414.89, -8419.06, -8362.64, -8381.55, -8370.53, -8392.66, -8389.77, -8393.96, -8388.11, -8381.73, -8333.95, -8348.93, -8363.93, -8360.30};
    double RS2_L_off[16]={-8473.92, -8470.25, -8489.89, -8494.06, -8512.64, -8494.05, -8520.53, -8505.16, -8502.27, -8468.96, -8500.61, -8494.23, -8521.45, -8498.93, -8476.43, -8472.80};
    double LS2_R_off[16]={-12441.14, -12490.70, -12579.43, -12601.39, -12471.56, -12471.38, -12658.08, -12656.28, -12690.65, -12489.77, -12701.56, -12675.30, -12696.36, -12491.35, -12709.36, -12539.99 };
    double LS2_L_off[16]={-12141.14, -12190.70, -12091.93, -12076.39, -12209.06, -12208.88, -12058.08, -12056.28, -12015.65, -12227.27, -12026.56, -12000.30, -11983.86, -12191.35, -11996.86, -12239.99 };

    double  RS2_off_H2[16];
    double  LS2_off_H2[16];
    for(int l=0;l<16;l++){
      RS2_off_H2[l]=RS2_R_off[l] + RS2_L_off[l];
      LS2_off_H2[l]=LS2_R_off[l] + LS2_L_off[l];
    }
      
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(KINE==1){

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,-16895.6,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}

