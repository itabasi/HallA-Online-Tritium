

void hrs_KID_F1_pi(){

  TChain*  T=new TChain("tree");
  tree->Add("/home/itabashi/jlab_nnL/rootfiles/Lambda_small_mm_ana1224_woAC.root");

  double ctime,sum_a1,sum_a2;
 T->SetBranchStatus("*",0);  
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&ctime);
 T->SetBranchStatus("ac1_sum",1);    
 T->SetBranchAddress("ac1_sum",&sum_a1);
 T->SetBranchStatus("ac2_sum",1);    
 T->SetBranchAddress("ac2_sum",&sum_a2);
 TH1F* hcoin_acc=new TH1F("hcoin_acc","Coin ACC region Hist",1000,-40,20.);
 TH1F* hacc=new TH1F("hacc","ACC Hist",1000,-40.,20.);
 TH1F* hcoin=new TH1F("hcoin","Coin time Hist",1000,-40,20.);
 TH1F* hcoin_pi=new TH1F("hcoin_pi","Coin time Hist",1000,-40,20.);

 int evnt=tree->GetEntries();
 cout<<"Event num :"<<evnt<<endl;
 // evnt=1000;
 double time;
 for(int k=0;k<evnt;k++){
   tree->GetEntry(k);
   hcoin->Fill(ctime);
   if(5000.<sum_a2 && 300<sum_a1){
   hcoin_pi->Fill(ctime);
   if((-63<ctime && ctime <-15) || (15<ctime && ctime<63)){
       double ct=ctime;       
        while(1){
	  if(-3.0<ct && ct<3.0){
		 hacc->Fill(ct);
                 hacc->Fill(ct-36);
		 break;}
	       else if(ct<-3.0){ct=ct+6;}
	       else if(3.0<ct){ct=ct-6;}
	 }
     }
   }
 }

     hacc->Scale(6.0/96.);

     TF1* fpi=new TF1("fpi","gausn(0)+gausn(3)",0.0,5.0);
     TF1* facc=new TF1("facc","gausn(0)",0.0,5.0);
     double pi[7],acc[7];
     pi[0]=4.60353e+04;
     pi[1]=3.05727e+00;
     pi[2]=2.61242e-01;
     pi[3]=1.46537e+03;     
     pi[4]=2.0;
     pi[5]=3.54211e-01;

     acc[0]=1.46537e+03;     
     acc[1]=1.0;
     acc[2]=3.54211e-01;

     fpi->SetParameters(pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);
     facc->SetParameters(acc[0],acc[1],acc[2]);
     facc->FixParameter(0,acc[0]);
     facc->FixParameter(1,acc[1]);
     facc->FixParameter(2,acc[2]);

     //  fpi->FixParameter(3,acc[0]);
     //     fpi->FixParameter(4,acc[1]);
     fpi->SetParLimits(4,acc[1]-0.5,acc[1]+0.5);
     //     fpi->FixParameter(5,acc[2]);
     hcoin_pi->Fit("fpi","","",0.0,4.0);

 TCanvas* c0=new TCanvas("c0","Coin time ACC region Hist");
 c0->cd();
 hcoin_pi->Draw();
 fpi->SetLineColor(2);
 fpi->SetFillColor(2);
 fpi->SetFillStyle(3001);
 fpi->Draw("same");
 facc->SetLineColor(3);
 facc->SetFillColor(3);
 facc->SetFillStyle(3001);
 facc->Draw("same");

 /*
 hacc->SetLineColor(4);
 hacc->SetFillColor(4);
 hacc->SetFillStyle(3001);
 hacc->Draw("same");
 */
 //TCanvas* c1=new TCanvas("c1","ACC Hist");
 // c1->cd();

      }
