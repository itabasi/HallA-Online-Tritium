/*
  singleout.cc
  Toshiyuki Gogami, Nov 22, 2019
*/


void singleout(){
  //ifstream* ifs = new ifstream("Rootfiles/h2.dat");ofstream* ofs = new ofstream("h2_replay.dat");
  //ifstream* ifs = new ifstream("Rootfiles/h22.dat");ofstream* ofs = new ofstream("h22_replay.dat");
  //ifstream* ifs = new ifstream("Rootfiles/T2.dat");ofstream* ofs = new ofstream("T2_replay.dat");
  ifstream* ifs = new ifstream("Rootfiles/He3.dat");ofstream* ofs = new ofstream("He3_replay.dat");
  
  int a, b, temp;
  *ifs >> a >> temp >> temp;
  b = a;
  *ofs << a << endl;
  while (!ifs->eof()){
    *ifs >> a >> temp >> temp;
    if(a!=b){
      cout << a << endl;
      *ofs << a << endl;
      b=a;
    }
  }
  ofs->close();
}
