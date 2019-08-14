
void test(){

  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  int  nMatT=5;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      if (a+b+c+d+e==n){
		cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl; 
		npar++;
	      }
	    }
	  }
	}
      }    
    }
  }


  cout<<"npar "<<npar<<endl;

}
