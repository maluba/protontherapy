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
// $Id: B4cEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "HEventAction.hh"
#include "phantomSD.hh"
#include "phantomHit.hh"
#include "HAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HEventAction::HEventAction()
 : G4UserEventAction(),
   fSliceZHCID(-1)
//   fCylHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HEventAction::~HEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomHitsCollection* 
HEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  phantomHitsCollection* hitsCollection 
    = static_cast<phantomHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("HEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HEventAction::PrintEventStatistics(
                              G4double sliceZEdep, G4double sliceZTrackLength
                             /*, G4double cylEdep, G4double cylTrackLength*/) const
{
  // print event statistics
//  G4cout
//     << "   SliceZ: total energy: " 
//     << std::setw(7) << G4BestUnit(sliceZEdep, "Energy")
//     << "       total track length: " 
//     << std::setw(7) << G4BestUnit(sliceZTrackLength, "Length")
//     << G4endl;

//     << "        SliceR: total energy: " 
//     << std::setw(7) << G4BestUnit(cylEdep, "Energy")
//     << "       total track length: " 
//     << std::setw(7) << G4BestUnit(cylTrackLength, "Length")
//     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
// G4int evtNb = evt->GetEventID();
// 
// //printing survey
// if (evtNb%fPrintModulo == 0) 
//    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HEventAction::EndOfEventAction(const G4Event* evt)
{  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  
  // Get hits collections IDs (only once)
  if ( fSliceZHCID == -1) {
    fSliceZHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSliceHitsCollection"); 
//    fCylHCID 
//      = G4SDManager::GetSDMpointer()->GetCollectionID("TrackerCylinderHitsCollection");
  }

  // Get hits collections
  phantomHitsCollection* sliceZHC = GetHitsCollection(fSliceZHCID, evt);

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//-----------------------------------------------------------------------
  G4double ZHalfLength = 0;
//  G4double RadiiHalfLength = 0;
  G4LogicalVolume* CylLV = G4LogicalVolumeStore::GetInstance()->GetVolume("RingLV");
  G4Tubs* slicethickness = NULL;
  if ( CylLV) slicethickness = dynamic_cast< G4Tubs*>(CylLV->GetSolid()); 
  if ( slicethickness ) {ZHalfLength = slicethickness->GetZHalfLength()/cm; } 
//  if ( slicethickness ) {RadiiHalfLength = slicethickness->GetXHalfLength()/cm; }

  else  {
    G4ExceptionDescription msg;
    msg << "volume not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
//    msg << "The gun will be place in the center.";
    G4Exception("HEventAction::EndOfEventAction()",
      "MyCode0002", JustWarning, msg);
  } 

  //Histogramming
  G4int nofCells = sliceZHC->entries();
  for (G4int i=0; i<nofCells; i++)
  {
   // if ((*sliceZHC)[i]->GetParticle()->GetParticleName() == "proton" && (*sliceZHC)[i]->GetTrackID() == 1 )
  // fill histgrams
//    {
//    G4cout << (*sliceZHC)[i]->GetPosition() << G4endl;
     // analysisManager->FillH1(1, (2*i+1)*ZHalfLength, (*sliceZHC)[i]->GetEdep()/MeV);
//    }
//    else if ((*sliceZHC)[i]->GetParticle()->GetParticleName() == "proton" && (*sliceZHC)[i]->GetTrackID() > 1 )
//    {
//      analysisManager->FillH1(2, (2*i+1)*ZHalfLength, (*sliceZHC)[i]->GetEdep()/MeV);
//    }
//    else {analysisManager->FillH1(3, (2*i+1)*ZHalfLength, (*sliceZHC)[i]->GetEdep()/MeV);}

    // fill histgrams
    analysisManager->FillH1(1, (2*i+1)*ZHalfLength, (*sliceZHC)[i]->GetEdep()/MeV);
    analysisManager->FillH1(2, (2*i+1)*ZHalfLength, (*sliceZHC)[i]->GetTrackLength()/cm);

    analysisManager->FillH1(3, (2*i+5)*(ZHalfLength), (*sliceZHC)[i]->GetEdep()/MeV);
    analysisManager->FillH1(4, (2*i+5)*(ZHalfLength), (*sliceZHC)[i]->GetEdep()/MeV);
//    analysisManager->FillH1(3, (*sliceZHC)[i]->GetEdep());
//    analysisManager->FillH1(3, (sliceZHC->GetEdep()));
  }  

 
  // Print per event (modulo n)
  //
    G4int eventID = evt->GetEventID();
    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
      G4cout << "---> End of event: " << eventID << G4endl;     

//      PrintEventStatistics(
//        (*sliceZHC)[i]->GetEdep(), (*sliceZHC)[i]->GetTrackLength());
//        (*cylHC)[i]->GetEdep(), (*cylHC)[i]->GetTrackLength());
  }  
  

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
