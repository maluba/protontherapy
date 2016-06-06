//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
/*
Author: Susanna Guatelli
*/
// The class BrachyAnalysisManager creates and manages histograms and ntuples

// The analysis was included in this application following the extended Geant4
// example analysis/AnaEx01

#include <stdlib.h>
#include "HAnalysis.hh"
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "time.h"

//#ifdef G4ANALYSIS_USE_ROOT
//#include "TROOT.h"
//#include "TFile.h"
//#include "TNtuple.h"
//#include "TH1F.h"
//#include "TH2F.h"
//#endif

AnalysisManager::AnalysisManager()
{

factoryOn = false;

// Initialization histograms
for (G4int k=0; k<MaxHisto; k++) {fHistId[k] = 0;}

//  // Creating histograms
//  analysisManager->CreateH1("1","Edep in absorber", 100, 0., 800*MeV);
//  analysisManager->CreateH1("2","Edep in gap", 100, 0., 100*MeV);
//  analysisManager->CreateH1("3","trackL in absorber", 100, 0., 1*m);
//  analysisManager->CreateH1("4","trackL in gap", 100, 0., 50*cm);


// Initialization ntuple
  for (G4int k=0; k<MaxNtCol; k++) 
  {
    fNtColId[k] = 0;
  }  

PrimaryProtonEdep = 0;
SecondaryProtonEdep = 0;
OtherSecondaryParticlesEdep = 0;

}

AnalysisManager::~AnalysisManager() 
{ 
}

void AnalysisManager::book() 
{  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  //enable activation of histograms
 
  // Create a root file
  G4String fileName = "hadrontherapy.root";

  // Create directories  
  analysisManager->SetHistoDirectoryName("hadrontherapy_histo");
  analysisManager->SetNtupleDirectoryName("hadrontherapy_ntuple");
  

  G4bool fileOpen = analysisManager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName[1] 
           << G4endl;
    return;
  }

//creating a 1D histograms
  analysisManager->SetFirstHistoId(1);
//  G4int fNofSlicesAlongZ = 300;
  
  // Histogram containing the primary proton energy (MeV) 
//  for (G4int k=0; k<MaxHisto; k++) { }
  fHistId[1]  = analysisManager -> CreateH1("1", "Initial energy", 1000, 0., 1000.);

  fHistId[2]  = analysisManager -> CreateH1("2", "Initial energy",  1000, 0., 1000.);

  fHistId[3]  = analysisManager -> CreateH1("3", "Initial energy", 1000, 0., 1000.);

  //Parameters of CreateH1: histoID, histo name, bins' number, xmin, xmax 
  
  PrimaryProtonEdep = analysisManager-> GetH1(fHistId[1]);
  SecondaryProtonEdep = analysisManager-> GetH1(fHistId[2]);
  OtherSecondaryParticlesEdep = analysisManager-> GetH1(fHistId[3]);
  
//  //creating a ntuple, containg 3D energy deposition in the phantom
//  analysisManager -> CreateNtuple("1", "3Dedep");
//  fNtColId[0] = analysisManager->CreateNtupleDColumn("xx");
//  fNtColId[1] = analysisManager->CreateNtupleDColumn("yy");
//  fNtColId[2] = analysisManager->CreateNtupleDColumn("zz");
//  fNtColId[3] = analysisManager->CreateNtupleDColumn("edep");
//  analysisManager->FinishNtuple();
  
 factoryOn = true;    
}

void AnalysisManager::FillPrimaryProtonHistogram(G4double primaryParticleEnergy)
{
//  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 // 1DHistogram: energy spectrum of primary particles  
  PrimaryProtonEdep -> fill(primaryParticleEnergy);
}
 
void AnalysisManager::FillScondaryProtonHistogram(G4double secondaryProtonEnergy)
{
 // 1DHistogram: energy spectrum of secondary protons  
  SecondaryProtonEdep -> fill(secondaryProtonEnergy);
}

void AnalysisManager::FillOtherSecondariesHistogram(G4double otherSecondariesEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  OtherSecondaryParticlesEdep -> fill(otherSecondariesEnergy);
}

//void AnalysisManager::FillNtupleWithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double zz,
//                                                     G4double energyDep)
//{
//  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//  analysisManager->FillNtupleDColumn(fNtColId[0], xx);
//  analysisManager->FillNtupleDColumn(fNtColId[1], yy);
//  analysisManager->FillNtupleDColumn(fNtColId[2], zz);
//  analysisManager->FillNtupleDColumn(fNtColId[3], energyDep);
//  analysisManager->AddNtupleRow();  
//}

void AnalysisManager::save() 
{  
 if (factoryOn) 
   {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
   }
}



