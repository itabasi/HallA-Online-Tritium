
void num_param_cal(){

  int nMatT=6;

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
			npar++;
		      }
		    }
		  }
		}
	      }    
	    }
  }

  cout<<"order "<<nMatT<<" npar "<<npar<<endl;

}

