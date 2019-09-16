// Aurther Itabshi 2019 June 30th
// LHRS Anguler resolution study in sieve slit data #111731
// ==========================


void angL_resolution(){


const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};

double selection_width = 0.01; // event selection width for z


const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider
 const double l0 = 100.3;
const double hrs_ang = 13.2 * 3.14159 / 180.;


  
  string buf;
  int hole, foil;
  double ssy_cent_real[nrow];
  double ssx_cent_real[ncol];
  double l[100],projectf[100],dth[100];

  for(int i=0;i<nfoil;i++){
    for(int j=0;j<nsshole;j++){

      
    l[i] = 0;
    l[i] = sqrt(pow(l0,2.0) + pow(fcent_real[i]*100.,2.0) -2.0*l0*fcent_real[i]*100.*cos(hrs_ang));
    dth[i] = asin(l0/l[i]*sin(hrs_ang)) - hrs_ang;
    projectf[i] = cos( dth[i] );
    cout<<"i "<<i<<" l "<<l[i]<<" dth "<<dth[i]*180./3.14<<" project "<<projectf[i]<<endl;
      }
  }
}
