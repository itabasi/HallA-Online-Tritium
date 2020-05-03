/*
  charge_int.cc
  
  Toshi Gogami, Mar 12, 2020
*/

void charge_int(){
  //ifstream* ifs = new ifstream("charge_T2.dat");
  //ifstream* ifs = new ifstream("charge_h2.dat");
  //ifstream* ifs = new ifstream("charge_h22.dat");
  ifstream* ifs = new ifstream("charge_He3.dat");
  int run;
  double charge;
  double sum = 0.0;
  while (!ifs->eof()){
    *ifs >> run >> charge;
    cout << " " << run << " " << charge << endl;
    sum = sum + charge;
  }
  cout << " ----------------" << endl;
  cout << " ==> Sum = " << sum/1.0e+6 << " C" << endl;
}

