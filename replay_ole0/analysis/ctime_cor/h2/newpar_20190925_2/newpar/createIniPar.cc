/*
  createIniPar.cc
  
  Toshiyuki Gogami, 25 Sep 2019
*/

void createIniPar(){
  
  const int nMatT=3;  
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  const int nZt=3;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  ofstream* ofs = new ofstream("ini.dat");
  
  double par=0;
  double temp=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      par=0;
	      temp = 0;
	      if (a+b+c+d+e==n){
		//if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){}
		//else{ par = 0.0;}
		cout << par << " " << a << " " << b << " " << c << " " << d << " " << e << endl;
		*ofs << par << " " << a << " " << b << " " << c << " " << d << " " << e << endl;
		npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }
  
  cout << endl;
  cout << npar << endl;
}


//////////////////////////////////////////////////
double calcf2t_4th(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt,
		   int flag=0)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		if(flag==2) Y += x*P[npar+126]; // LHRS
		else Y += x*P[npar];    // RHRS
	      npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }
  
  return Y; 
  
}

