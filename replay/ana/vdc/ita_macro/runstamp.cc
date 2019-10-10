//////////////////////////////////////////////
//Auther Itabashi 2019/10/03
//Write Run time stamp for parameters optics
/////////////////////////////////////////////

#include<string>
#include<stdio.h>

void runstamp(){

  string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/DB/db_run.dat";
  ifstream ifrun(fname.c_str());  
  string buf;
  int nmax=5000;
  string date[nmax];
  string run[nmax];
  int iii=0;
  string p;
  if (ifrun.fail()){ cerr << "failed open files" <<fname.c_str()<<endl; exit(1);}

  
    while( getline(ifrun,buf) ){
      if(buf[0]=='-'){
	date[iii]=buf;
	ifrun >> p >> p >> p >> p >> p >>run[iii]; 
      
	//      cout<<"date "<<date[iii]<<endl;
	//	cout<<"run "<<run[iii]<<endl;
      }

      
      if(atoi(run[iii].c_str())>111111)      iii++;
      if(iii>nmax)break;
    }


    ifrun.close();

    
    int run_init,run_end;
    run_init=0; run_end=0;
    char c[100];
    string st;
    string init,end;
    //    string paraname="./param/111500-111509.dat";
    string paraname="./param/initial_10run/111500-111509.dat";
    paraname="./param/initial_10run/111510-111519.dat";
    paraname="./param/initial_all/VDCt0_all.dat";
    bool init_flag=false;
    int run_test;
    for(int i=0;i<(int)paraname.length();i++){
      st=paraname.substr(i,1);

      if(st=="1"){
	if(init_flag==0){
	  init=paraname.substr(i,6);
	  run_test=atoi(init.c_str());
	  if(run_test>111111){
	    run_init=run_test;
	    init_flag=true;
	  }
	}else{
	  init=paraname.substr(i,6);
	  run_test=atoi(init.c_str());
	  if(run_test>111111)	run_end=run_test;
	}
      }
    }

    cout<<"run_int "<<run_init<<" run_end "<<run_end<<endl;

    bool all_flag=false;
    all_flag=true;
    if(all_flag){
      run_init=111140;
      run_end=111840;
    }

    string ofname;
    bool RVDC_flag=true;
    RVDC_flag=false;    
    if(RVDC_flag)ofname="./DB/db_R.vdc.dat";
    else ofname="./DB/db_L.vdc.dat";
    ofstream ofs(ofname.c_str(),std::ios::app);
   
    ifstream ifp(paraname.c_str());
    if (ifp.fail()){ cerr << "failed open files" <<paraname<<endl; exit(1);}    
    int plane=-1;
    // plane 0 : RVDC U1  plane 5 : LVDC U1
    // plane 1 : RVDC U2  plane 6 : LVDC U2 
    // plane 2 : RVDC V1  plane 7 : LVDC V1
    // plane 3 : RVDV V2  plane 8 : LVDC V2
    
    int ii=0;
    const int nwire=368;
    double off[nwire][8];   
    while( getline(ifp,buf) ){
      
      if( buf[0]=='#' ){
	ii=0;
	plane++;}
      if( ifp.eof() || ii+1 > nwire){continue;}
      
      ifp >> off[ii][plane] >> off[ii+1][plane] >> off[ii+2][plane] >> off[ii+3][plane]
	      >> off[ii+4][plane] >> off[ii+5][plane] >> off[ii+6][plane] >> off[ii+7][plane];
      ii=ii+8;
    }//end while dat file


    //==== Change plane ======//

    double offset[nwire][8];
    // plane 0 : RVDC U1  plane 5 : LVDC U1
    // plane 1 : RVDC V1  plane 6 : LVDC V1
    // plane 2 : RVDC U2  plane 7 : LVDC U2     
    // plane 3 : RVDV V2  plane 8 : LVDC V2

    for(int i=0;i<nwire;i++){
      offset[i][0]=off[i][0]; //RVDC U1
      offset[i][1]=off[i][2]; //RVDC V1
      offset[i][2]=off[i][1]; //RVDC U2
      offset[i][3]=off[i][3]; //RVDC V2
      offset[i][4]=off[i][4]; //LVDC U1
      offset[i][5]=off[i][6]; //LVDC V1
      offset[i][6]=off[i][5]; //LVDC U2
      offset[i][7]=off[i][7]; //LVDC V2      
      
      

    } 
    
    //==== Get date of run_init ===//
    string vdc[8];
    vdc[0]="R.vdc.u1.tdc.offsets = "; // RVDC U1
    vdc[1]="R.vdc.v1.tdc.offsets = "; // RVDC V1
    vdc[2]="R.vdc.u2.tdc.offsets = "; // RVDC U2
    vdc[3]="R.vdc.v2.tdc.offsets = "; // RVDC V2
    vdc[4]="L.vdc.u1.tdc.offsets = "; // LVDC U1
    vdc[5]="L.vdc.v1.tdc.offsets = "; // LVDC V1
    vdc[6]="L.vdc.u2.tdc.offsets = "; // LVDC U2
    vdc[7]="L.vdc.v2.tdc.offsets = "; // LVDC V2        


    string Date;
    
    for(int i=0;i<nmax;i++){
      if(atoi(run[i].c_str())==run_init){
	Date=date[i];	break;}
    }
    
    ofs<<Date.c_str()<<endl;
    for(int j=0;j<4;j++){
      if(RVDC_flag)ofs<<vdc[j]<<endl;
      else ofs<<vdc[j+4]<<endl;
      int trans=0;
      for(int i=0;i<nwire/8;i++){
	
	if(RVDC_flag)	ofs <<offset[i][j]<<" "<<offset[i+1][j]<<" "<<offset[i+2][j]<<" "<<offset[i+3][j]<<" "<<
			  offset[i+4][j]<<" "<<offset[i+5][j]<<" "<<offset[i+6][j]<<" "<<offset[i+7][j]<<" "<<endl;

	else 	ofs <<offset[i][j+4]<<" "<<offset[i+1][j+4]<<" "<<offset[i+2][j+4]<<" "<<offset[i+3][j+4]<<" "<<
			  offset[i+4][j+4]<<" "<<offset[i+5][j+4]<<" "<<offset[i+6][j+4]<<" "<<offset[i+7][j+4]<<" "<<endl;
      }
      ofs<<endl;
    }

    ifp.close();
    ofs.close();
    
    
}
