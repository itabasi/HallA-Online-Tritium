//////////////////////////
// AC gain parameter optimization macro
// Auther itabashi Nov. 19th ,2019
//////////////////////////

extern double gain2(int seg);
extern double gain3(int seg);
void ac_gain(){

  //  TChain* oldtree=new TChain("T");
  int ENum;
  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH1_1113.root";
  TChain *T = new TChain("T");
  T->Add(ifname.c_str());

  int nac1=24;
  int nac2=26;
  double ac1_npe[nac1],ac2_npe[nac2],ac1_sum,ac2_sum;
  double ac1_npe_sum,ac2_npe_sum,mm,ct,Rx_fp;
  int z_cut,ct_cut, Rs2_pad[100],Ls2_pad[100];
  //  double ac1_npe_c[nac1],ac2_npe_c[nac2];
  
  T->SetBranchStatus("*",0);
  // ana_Lambda root //
  T->SetBranchStatus("mm",1);
  T->SetBranchAddress("mm",&mm);
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);  
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);
  T->SetBranchStatus("ct_cut",1);
  T->SetBranchAddress("ct_cut",&ct_cut);      
  T->SetBranchStatus("ct",1);
  T->SetBranchAddress("ct",&ct);
  T->SetBranchStatus("Rs2_pad",1);
  T->SetBranchAddress("Rs2_pad",Rs2_pad);
  T->SetBranchStatus("Ls2_pad",1);
  T->SetBranchAddress("Ls2_pad",Ls2_pad);
  
  T->SetBranchStatus("ac1_sum",1);
  T->SetBranchAddress("ac1_sum",&ac1_sum);
  T->SetBranchStatus("ac2_sum",1);
  T->SetBranchAddress("ac2_sum",&ac2_sum);
  T->SetBranchStatus("ac1_sum",1);
  T->SetBranchAddress("ac1_sum",&ac1_sum);
  T->SetBranchStatus("ac2_sum",1);
  T->SetBranchAddress("ac2_sum",&ac2_sum);
  T->SetBranchStatus("ac1_npe",1);
  T->SetBranchAddress("ac1_npe",ac1_npe);
  T->SetBranchStatus("ac2_npe",1);
  T->SetBranchAddress("ac2_npe",ac2_npe);
  T->SetBranchStatus("ac2_npe",1);
  T->SetBranchAddress("ac2_npe",ac2_npe);  
  T->SetBranchStatus("ac1_npe_sum",1);
  T->SetBranchAddress("ac1_npe_sum",&ac1_npe_sum);
  T->SetBranchStatus("ac2_npe_sum",1);
  T->SetBranchAddress("ac2_npe_sum",&ac2_npe_sum);  

  string ofname="../rootfiles/ac/test.root";
  TFile* ofp= new TFile(ofname.c_str(),"recreate");
  TTree* tnew=new TTree("T","AC gain optimization");
  tnew=T->CloneTree(0);

  double ac1_npe_c[nac1],ac2_npe_c[nac2],ac1_npe_sum_c,ac2_npe_sum_c;
  tnew->Branch("ac1_npe_c",ac1_npe_c,"ac1_npe_c[24]/D");
  tnew->Branch("ac2_npe_c",ac2_npe_c,"ac2_npe_c[26]/D");
  tnew->Branch("ac1_npe_sum_c",&ac1_npe_sum_c,"ac1_npe_sum_c/D");  
  tnew->Branch("ac2_npe_sum_c",&ac2_npe_sum_c,"ac2_npe_sum_c/D");  

  ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

  for(int k=0;k<ENum;k++){
    T->GetEntry(k);

    ac2_npe_sum_c=0.0;
    for(int i=0;i<nac1;i++)ac1_npe_c[i]=-100.;
        for(int i=0;i<nac2;i++)ac2_npe_c[i]=-100.;

	/*
	for(int i=0;i<nac2;i++){
	  ac2_npe_c[i]=ac2_npe[i]/gain2(i);
	  ac2_npe_sum_c+=ac2_npe_c[i];
	  }
*/
	

	  ac2_npe_c[Rs2_pad[0]]=ac2_npe_sum/gain3(Rs2_pad[0])*20.0;
	  ac2_npe_sum_c+=ac2_npe_c[Rs2_pad[0]];


	tnew->Fill();
  }

}


double gain2(int seg){

  double gain[26]={
		   1.22680e+00,
		   9.54319e-01,
		   9.58565e-01,
		   3.17281e-01,
		   9.69098e-01,
		   5.24055e-01,
		   9.47421e-01,
		   9.54565e-01,
		   9.75232e-01,
		   9.69701e-01,
		   9.67153e-01,
		   1.00151e+00,
		   9.83539e-01,
		   9.86966e-01,
		   9.91463e-01,
		   1.00155e+00,
		   9.81375e-01,
		   1.00786e+00,
		   9.24110e-01,
		   9.68688e-01,
		   9.41630e-01,
		   9.67952e-01,
		   1.02218e+00,
		   9.93827e-01,
		   1.06834e+00,
		   9.81790e-01,

		                       

  };

    return gain[seg];
  

}

double gain3(int seg){

  double gain[16]={1.0,
		   2.07850e+01,
		   1.92837e+01,
		   1.83174e+01,
		   1.76284e+01,
		   1.80764e+01,
		   1.75492e+01,
		   1.96382e+01,
		   2.01166e+01,
		   1.95302e+01,
		   2.14401e+01,
		   2.03358e+01,
		   1.82856e+01,
		   1.62902e+01,
		   1.43125e+01,
		   1.34920e+01

  };

  return gain[seg];
}
