// GetRootFileName.C
//
// Helper macro for OnlineGUI to determine ROOT Filename
//   from a given runnumber
//
// 
TString GetRootFileName(UInt_t runnumber)
{
  TString file_path = "ROOTfiles/pass1";
  TString expID = "tritium_online";
  
  TString filename;
  filename = file_path + expID;
  filename += "_";
  filename += runnumber;
  filename += "-all.root";

  return filename;

}
