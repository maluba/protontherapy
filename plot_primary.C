{
// Plot the energy spectrum of the primary particles
gROOT -> Reset();

TFile f("hadrontherapy.root");
  
// Draw histos filled by Geant4 simulation 
//   
  
TCanvas* c1 = new TCanvas("c1", " ");

TDirectory* dir_histo = f.Get("hadrontherapy_histo");
    

// Primary particle spectrum
TH1D* hist1 = (TH1D*)dir_histo->Get("1");   
hist1->Draw(" "); // this is empy at the moment.
}
