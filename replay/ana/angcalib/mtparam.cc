//////// function

double comb(int n, int k){

  double y;
  
  int x0=1;
  int x1=1;

  if(k==0)return 1;
  for(int i=n-k+1 ; i<=n; i++ )x0 = x0*i;
  for(int i=1;      i<=k; i++ )x1 = x1*i;
  y=(double)x0/(double)x1;
  return y;

}




void mtparam(){

  //  string fname="./db_RHRS_theta.dat";
  //  string ofname="./test_xp.dat";
  //  string fname="./db_RHRS_phi.dat";
  //  string ofname="./test_yp.dat";
  string fname="./db_RHRS_y.dat";
  string ofname="./test_y.dat";  
  ofstream * ofs1=new ofstream(ofname.c_str());

  ifstream mtp(fname.c_str());
  string buf;

  string type;
  int a,b,c;
  double x0,x1,x2,x3,x4;
  int nn=5;
  double par[nn][nn][nn][nn][nn];
  int i=0;


const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  PLm = 25.4, PLr=0.7; 
const double  Ztm = -0.15,Ztr=0.35; 

  
  for(int i=0 ; i<nn ; i++){
      for(int e=0 ; e<nn ; e++){
	for(int d=0 ; d<nn ; d++){
	  for(int c=0 ; c<nn ; c++){
	    for(int b=0 ; b<nn ; b++){
	      for(int a=0 ; a<nn ; a++){  
		par[a][b][c][d][e]=0.0;
	      }
	    }
	  }
	}
      }
  }

  double x[5],y[5],th[5],ph[5];
  while( getline(mtp,buf) ){
    a=0; b=0; c=0; x0=0.0; x1=0.0; x2=0.0; x3=0.0; x4=0.0;
    sscanf(buf.c_str(),"  %s %d %d %d %lf %lf %lf %lf %lf",type.c_str(),&a, &b, &c, &x0, &x1, &x2, &x3, &x4);
    cout<<type.c_str()<<" "<<a<<" "<<b<<" "<<c<<" "<<x0<<" "<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<endl;


    par[a][b][c][0][0]=x0;
    par[a][b][c][1][0]=x1;
    par[a][b][c][2][0]=x2;
    par[a][b][c][3][0]=x3;
    par[a][b][c][4][0]=x4;

    
    /*
    par[0][a][b][c][0]=x0;
    par[1][a][b][c][0]=x1;
    par[2][a][b][c][0]=x2;
    par[3][a][b][c][0]=x3;
    par[4][a][b][c][0]=x4;
    */

  }// end while


    int nppp = 0;

    for(int i=0 ; i<nn ; i++){
      for(int e=0 ; e<nn ; e++){
	for(int d=0 ; d<nn ; d++){
	  for(int c=0 ; c<nn ; c++){
	    for(int b=0 ; b<nn ; b++){
	      for(int a=0 ; a<nn ; a++){  
		if(a+b+c+d+e==i){
		  *ofs1 << par[a][b][c][d][e]
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  //		  cout<<par[a][b][c][d][e]<<endl;
		}
	      }
	    }
	  }
	}
      }
    }
     

    ofs1->close();
    
}


