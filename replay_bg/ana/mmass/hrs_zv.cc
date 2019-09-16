
void hrs_zv(){

  char* hyd_root="../../rootfiles/Lambda_small1_mm_ana1224_woAC.root";
  TFile *fin_hyd =new TFile(hyd_root,"read");

  char* emp_root="../../rootfiles/tritium_111325.root";
  TFile *fin_emp =new TFile(emp_root,"read");

  TTree *t_hyd=(TTree*)fin_hyd->Get("tree");
  TTree *t_emp=(TTree*)fin_emp->Get("T");
  double Rz_hyd,Rz_emp[100];
  t_hyd->SetBranchAddress("Rz",&Rz_hyd);
  t_emp->SetBranchStatus("*",0);  
  t_emp->SetBranchStatus("R.tr.vz",1);
  t_emp->SetBranchAddress("R.tr.vz",Rz_emp);
  TH1F* hvz_hyd=new TH1F("hvz_hyd","Hydrogen run",1000,-0.2,0.2);
  TH1F* hvz_emp=new TH1F("hvz_emp","Empty run",1000,-0.2,0.2);

  int ev_hyd=t_hyd->GetEntries();
  int ev_emp=t_emp->GetEntries();

  //  ev_hyd=100000.;
  // ev_emp=10000;
 int k=0;
  for(int i=0;i<ev_hyd;i++){
    t_hyd->GetEntry(i);
    if(i==10000*k){cout<<"i :"<<i<<"/"<<ev_hyd<<endl; k=k+1; }
    hvz_hyd->Fill(Rz_hyd);}
 
  k=0;
  for(int i=0;i<ev_emp;i++){
    t_emp->GetEntry(i);
    if(i==10000*k){cout<<"i :"<<i<<"/"<<ev_emp<<endl; k=k+1; }
    hvz_emp->Fill(Rz_emp[0]);}



  cout<<"Fill is done !!"<<endl;

    TF1* fhyd=new TF1("fhyd","gaus(0)",0.1,0.2);
    TF1* femp=new TF1("femp","gaus(0)",0.1,0.2);

    double ph[3],pe[3];
    ph[0]=1.35463e+03;
    ph[1]=1.37262e-01;
    ph[2]=8.29885e-03;
    fhyd->SetParameters(ph[0],ph[1],ph[2]);
    pe[0]=1.35463e+03;
    pe[1]=1.37262e-01;
    pe[2]=8.29885e-03;
    femp->SetParameters(pe[0],pe[1],pe[2]);
   hvz_hyd->Fit("fhyd","RqN","",0.13,0.15);
   hvz_emp->Fit("femp","RqN","",0.13,0.15);

  cout<<"Fitting is done !!"<<endl;


  double p0_hyd=fhyd->GetParameter(0);
  double p0_emp=femp->GetParameter(0);


  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  
  hvz_hyd->SetLineColor(2);
  hvz_hyd->SetFillColor(2);
  hvz_hyd->SetFillStyle(3001);
  hvz_emp->SetLineColor(4);
  hvz_emp->SetFillColor(4);
  hvz_emp->SetFillStyle(3001);

  hvz_hyd->Scale(p0_emp/p0_hyd);
  hvz_hyd->Draw();
  hvz_emp->Draw("same");

TCanvas* c1 =new TCanvas("c1","c1");
 c1->Divide(1,2);
 c1->cd(1);
  hvz_hyd->Draw();
 c1->cd(2);
  hvz_emp->Draw();


}
