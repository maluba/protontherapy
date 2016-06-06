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
//
// $Id: B4RunAction.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "HRunAction.hh"
#include "HAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HRunAction::HRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
//  G4double SliceThicknessZ = 2 *mm;
  G4int fNofSlices = 60;
  G4int fNbofRings = 10;
  
  // Creating histograms
//  analysisManager->CreateH1("1","Primary Proton Edep", fNofSlices, 0., 30);
//  analysisManager->CreateH1("2","Secondary Proton Edep", fNofSlices, 0., 30);
//  analysisManager->CreateH1("3","Recoil ions and other Secondaries Edep", fNofSlices, 0., 30);
//  analysisManager->CreateH1("1","Edep", fNofSlices, 0., 30);
  
  analysisManager->CreateH1("1","Edep", fNofSlices, 0., 30);
  analysisManager->CreateH1("2","trackL", fNofSlices, 0., 30);
 
  analysisManager->CreateH1("3","Lateral Edep", fNbofRings, 0., 30);
  analysisManager->CreateH1("4","trackLength in radial distance", fNbofRings, 0., 30);
  
//  analysisManager->CreateH1("3","dose frequency along Z-axis", fNofSlicesAlongZ, 0., 20*MeV);
//  analysisManager->CreateH1("4","tracklenghth frequency along Z-axis", NofSlicesAlongZ, 0., 30*cm);
  
//    analysisManager->CreateH1("3","Edep in cylinder along the z-axis", 100, 0., 100*MeV);
  // Creating ntuple
  //
/*  analysisManager->CreateNtuple("hadrontherapy", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("EsliceZ");
  analysisManager->CreateNtupleDColumn("EsliceR");
  analysisManager->CreateNtupleDColumn("LslicesZ");
  analysisManager->CreateNtupleDColumn("LslicesR");
  analysisManager->FinishNtuple();
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HRunAction::~HRunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "hadrontherapy";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if ( nofEvents == 0 ) return;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << "---> normalizing Hist  "<< G4endl; 
  G4double fac = (1./(nofEvents));
  if ( analysisManager->GetH1(1) ) analysisManager->GetH1(1)->scale(fac);
  if ( analysisManager->GetH1(2) ) analysisManager->GetH1(2)->scale(fac);
  if ( analysisManager->GetH1(3) ) analysisManager->GetH1(3)->scale(fac);
  if ( analysisManager->GetH1(4) ) analysisManager->GetH1(4)->scale(fac);

  
  // print histogram statistics
  //
//  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run \n" << G4endl; 
    }
    else {
      G4cout << "for the local thread \n" << G4endl; 
    }
    
//    G4cout << " ESliceZ : mean = " 
//       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
//       << " rms = " 
//       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    /*
    G4cout << " ESliceR : mean = " 
       << G4BestUnit(analysisManager->GetH1(3)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Energy") << G4endl;
    */
//    G4cout << " LSliceZ : mean = " 
//      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length") 
//      << " rms = " 
//      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;
/*
    G4cout << " LSliceR : mean = " 
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
*/
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
