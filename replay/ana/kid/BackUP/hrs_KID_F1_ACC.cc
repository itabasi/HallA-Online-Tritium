

void hrs_KID_F1_ACC(){

  TChain*  T=new TChain("tree");
  tree->Add("/home/itabashi/jlab_nnL/ita_scripts/rootfiles/coin_H2_1.root");

  double ctime[100],sum_a1,sum_a2;
 T->SetBranchStatus("*",0);  
 T->SetBranchStatus("ctime",1);    
 T->SetBranchAddress("ctime",ctime);
 T->SetBranchStatus("R.a1.asum_c",1);    
 T->SetBranchAddress("R.a1.asum_c",&sum_a1);
 T->SetBranchStatus("R.a2.asum_c",1);    
 T->SetBranchAddress("R.a2.asum_c",&sum_a2);
 TH1F* hcoin_acc=new TH1F("hcoin_acc","Coin ACC region Hist",1000,-20,20.);
 TH1F* hacc=new TH1F("hacc","ACC Hist",1000,-20.,20.);
 TH1F* hcoin=new TH1F("hcoin","Coin time Hist",1000,-20,20.);


 int evnt=tree->GetEntries();
 cout<<"Event num :"<<evnt<<endl;
 // evnt=1000;
 double time;

 for(int k=0;k<evnt;k++){
   tree->GetEntry(k);
   time=ctime[0];
   hcoin->Fill(time);

   if((-19< time && time < -4.) || 14.< time){
  
      while(1){
	   if(-19.<time && time<-17.){
	     hacc->Fill(time);
	     break;
	   }else{time=time-2.0; }
	 }
       }
 
   //   if(2.0<sum_a2 && sum_a2<7.0)hcoin_acc->Fill(ctime[0]);
       
 }

 int max;
 max=3;
 TF1* facc[1000];

 for(int i=0;i<max;i++){
 facc[i]=new TF1(Form("facc[%d]",i),"[0]+gaus(1)",-19+2*i,-17+2*i);
 facc[i]->SetParameter(0,60.);
 facc[i]->SetParLimits(1,10,100);
 facc[i]->SetParameter(2,-16.+2*i);
 facc[i]->SetParLimits(3,0.1,0.6);
 hcoin->Fit(Form("facc[%d]",i),"","",-19+2*i,-17+2*i);
 }


 TCanvas* c0=new TCanvas("c0","Coin time ACC region Hist");
 c0->cd();
 hcoin->Draw();
 for(int i=0;i<max;i++){
 facc[i]->SetLineColor(i+1);
 facc[i]->Draw("same");
 }

 TCanvas* c1=new TCanvas("c1","ACC Hist");
 c1->cd();
 hcoin->Draw();
 hacc->Scale(2./21);
 hacc->SetFillStyle(3001);
 hacc->SetFillColor(4);
 hacc->Draw("same");

      }
