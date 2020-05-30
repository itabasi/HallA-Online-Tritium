

void short_z(){


  string ifname ="../rootfiles/zcalib/multi_foil_old.root";
  TChain* t1=new TChain("T");
  t1->Add(ifname.c_str());  
  t1->SetBranchStatus("*",0);  
  t1->SetBranchStatus("fEvtHdr.fRun",1);
  t1->SetBranchStatus("HALLA_p",     1);
  t1->SetBranchStatus("DR.T1",       1);
  t1->SetBranchStatus("DR.T4",       1);
  t1->SetBranchStatus("DR.T5",       1);
  t1->SetBranchStatus("Ndata.R.tr.vz", 1);
  t1->SetBranchStatus("Ndata.L.tr.vz", 1);
  t1->SetBranchStatus("R.tr.vz",     1);
  t1->SetBranchStatus("L.tr.vz",     1);
  t1->SetBranchStatus("R.tr.x",      1);
  t1->SetBranchStatus("L.tr.x",      1);
  t1->SetBranchStatus("R.tr.y",      1);
  t1->SetBranchStatus("L.tr.y",      1);
  t1->SetBranchStatus("R.tr.vx",     1);
  t1->SetBranchStatus("L.tr.vx",     1);
  t1->SetBranchStatus("R.tr.vy",     1);
  t1->SetBranchStatus("L.tr.vy",     1);
  t1->SetBranchStatus("R.tr.th",     1);
  t1->SetBranchStatus("L.tr.th",     1);
  t1->SetBranchStatus("R.tr.ph",     1);
  t1->SetBranchStatus("L.tr.ph",     1);
  t1->SetBranchStatus("R.tr.tg_th",  1);
  t1->SetBranchStatus("R.tr.tg_ph",  1);
  t1->SetBranchStatus("L.tr.tg_th",  1);
  t1->SetBranchStatus("L.tr.tg_ph",  1);
  t1->SetBranchStatus("R.ps.asum_c", 1);
  t1->SetBranchStatus("R.sh.asum_c", 1);
  t1->SetBranchStatus("L.cer.asum_c",1);
  t1->SetBranchStatus("R.a1.asum_c", 1);
  t1->SetBranchStatus("R.a2.asum_c", 1);
  t1->SetBranchStatus("Lrb.Raster2.rawcur.x", 1); // raster current
  t1->SetBranchStatus("Lrb.Raster2.rawcur.y", 1); // raster current
  t1->SetBranchStatus("Rrb.Raster2.rawcur.x", 1); // raster current
  t1->SetBranchStatus("Rrb.Raster2.rawcur.y", 1); // raster current
  /*
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum);
  t1->SetBranchAddress("HALLA_p", &hallap);
  t1->SetBranchAddress("DR.T1", &trig1);
  t1->SetBranchAddress("DR.T4", &trig4);
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("Ndata.R.tr.vz",&Nrvz);
  t1->SetBranchAddress("Ndata.L.tr.vz",&Nlvz);
  t1->SetBranchAddress("R.tr.vz", rvz);
  t1->SetBranchAddress("L.tr.vz", lvz);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.vx",   &rx);
  t1->SetBranchAddress("L.tr.vx",   &lx);
  t1->SetBranchAddress("R.tr.vy",   &ry);
  t1->SetBranchAddress("L.tr.vy",   &ly);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.tr.tg_th", &Rth);
  t1->SetBranchAddress("R.tr.tg_ph", &Rph);
  t1->SetBranchAddress("L.tr.tg_th", &Lth);
  t1->SetBranchAddress("L.tr.tg_ph", &Lph);
  t1->SetBranchAddress("R.ps.asum_c", &R_ps_asum);
  t1->SetBranchAddress("R.sh.asum_c", &R_sh_asum);
  t1->SetBranchAddress("L.cer.asum_c",&L_gs_asum);
  t1->SetBranchAddress("R.a1.asum_c", &R_a1_asum);
  t1->SetBranchAddress("R.a2.asum_c", &R_a2_asum);
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &L_Ras_x); // raster current
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &L_Ras_y); // raster current
  t1->SetBranchAddress("Rrb.Raster2.rawcur.x", &R_Ras_x); // raster current
  t1->SetBranchAddress("Rrb.Raster2.rawcur.y", &R_Ras_y); // raster current  
  */
  string ofname ="./multi_foil_woras_small_old.root";                              
  TFile* ofs=new TFile(ofname.c_str(),"recreate");                             
  TTree* tnew = new TTree("T","Multi foilrun short root");                                     
  tnew =t1->CloneTree();
  tnew->Write();
}
