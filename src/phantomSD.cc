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
// $Id: B4cCalorimeterSD.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file B4cCalorimeterSD.cc
/// \brief Implementation of the B4cCalorimeterSD class

#include "phantomSD.hh"
#include "phantomHit.hh"
#include "HRunAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

  using namespace std;
//  ofstream myfile("myfilename");
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomSD::phantomSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomSD::~phantomSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void phantomSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new phantomHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  if(hcID<0)
    { hcID = GetCollectionID(0); }
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells more for total sums
  for (G4int i=0; i<fNofCells; i++ ) {
    fHitsCollection->insert(new phantomHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool phantomSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
//_________________________________________________________________________________________________________________
//  if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "PhantomSliceHitsCollection") return false;
//  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

//    // Get kinetic energy
    G4Track * theTrack = aStep  ->  GetTrack();
    G4double kineticEnergy =  theTrack -> GetKineticEnergy()/MeV;  

    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
//    //Get particle name  
    G4String particleName =  particleDef -> GetParticleName();  

    G4int pdg = particleDef ->GetPDGEncoding();

    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();

    // Energy Deposited
    G4double edep  = aStep -> GetTotalEnergyDeposit()/MeV;  //energyDeposit

    G4double stepLength  = aStep -> GetStepLength()/cm;  //DX

//    G4int Z = particleDef-> GetAtomicNumber();
//    G4int A = particleDef-> GetAtomicMass();
//  
    if ( edep==0. && stepLength == 0. ) return false;
  
    G4int pcharge = particleDef->GetPDGCharge();  //aStep->GetTrack()->GetDefinition()
    G4ThreeVector position = theTrack->GetPosition()/cm; //aStep->GetTrack()
    G4ThreeVector momentum = aStep->GetTrack()->GetMomentum()/MeV;
    G4String particleType = particleDef->GetParticleType();

  // Data collection
  //
  G4double sp,pp,rc;

  G4double x = position.x()/cm;
  G4double y = position.y()/cm;
  G4double z = position.z()/cm;

  G4double r = std::sqrt(x*x + y*y);
/*
  std::ofstream myfile ("tt.txt", ios::out | ios::app);
  if (myfile.is_open())
  {
   myfile << std::left << std::setw(20) << r << std::setw(20) << z << std::setw(20) << edep << G4endl;

   myfile.close();
  }
 */ 
  //---------------------------Edep due to particle category------------------------------

 if(particleName == "proton")
  {
    if ( aStep->GetTrack()->GetParentID() == 0 ) //applicable only to prim.
    {
      pp = aStep -> GetTotalEnergyDeposit()/MeV;
/*      std::ofstream myfile ("pp.txt", ios::out | ios::app);
      if (myfile.is_open())
    {
       myfile << std::left << std::setw(20) << r << std::setw(20) << z << std::setw(20) << edep << G4endl;

       myfile.close();
      }*/
      
    }
    else if( aStep->GetTrack()->GetParentID() != 0 ) // This is true only for sec.
    {
     sp = aStep -> GetTotalEnergyDeposit()/MeV;
/*      std::ofstream myfile ("sp.txt", ios::out | ios::app);
      if (myfile.is_open())
      {
       myfile << std::left << std::setw(20) << r << std::setw(20) << z << std::setw(20) << edep << G4endl;

       myfile.close();
      }*/
    }

//    if(PrmprotEdep == SecprotEdep)
//    {
//      SecprotEdep += edep;
//    }
  }
  else if(particleName != "proton")
  {
   rc = aStep -> GetTotalEnergyDeposit()/MeV;
/*      std::ofstream myfile ("rc.txt", ios::out | ios::app);
      if (myfile.is_open())
      {
       myfile << std::left << std::setw(20) << r << std::setw(20) << z << std::setw(20) << edep << G4endl;

       myfile.close();
      }*/
  }
  

//------------------------Secondaries-------------------------------------------

//          phantomHit* newHit = new phantomHit();
//          newHit->SetStripNo(  touchable->GetReplicaNumber(0) );
//          newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
//          newHit->SetMomentum( aStep->GetPreStepPoint()->GetMomentum() );
//          newHit->SetEdep( aStep->GetPreStepPoint()->GetTotalEnergy() );
//          newHit->SetTrackLength( aStep->GetStepLength() );                                    //GetStepLength() );
//          newHit->SetParticle( aStep->GetTrack()->GetDefinition() );
//          fHitsCollection->insert( newHit );

//_________________________________________________________________________________________________________________


//  G4String volumeName = aStep -> GetPreStepPoint() -> GetPhysicalVolume()-> GetName();
//  if(volumeName != "RingPhys")
//    return false;

//  phantomHit* newHit = new phantomHit();
//  newHit ->SetEdep(edep);
//  fHitsCollection ->insert(newHit);
//    //above is my new approach
    
  G4int sliceNumber = touchable->GetCopyNumber(0);
//  G4int ringNumber = touchable->GetCopyNumber(1);
  
//   Get hit accounting data for this cell
  phantomHit* hit = (*fHitsCollection)[sliceNumber];
 
  if ( !hit) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << sliceNumber; 
    G4Exception("phantomSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         
   

//   Add values
  hit->Add(sp, stepLength);

//  hit1->Add(PrmprotEdep, PrmTrkl);
//  hit2->Add(SecprotEdep, SecPartTrkl);
//  hit3->Add(OtherpartEdep, OtherpartTrkl);
   
   
  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void phantomSD::EndOfEvent(G4HCofThisEvent*)
//{
//  G4double totalEdepInOneEvent = 0;

//  if ( verboseLevel>1 ) { 
//     G4int nofHits = fHitsCollection->entries();
//     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
//            << " hits in the tracker slices: " << G4endl;
//     for ( G4int i=0; i<nofHits; i++ ) //(*fHitsCollection)[i]->Print();
//       {
//       G4double edep = (*fHitsCollection)[i] -> GetEdep();
//       totalEdepInOneEvent = totalEdepInOneEvent + edep;
//       }
//  if(totalEdepInOneEvent != 0) G4cout << "Total Energy per Event" << totalEdepInOneEvent << G4endl;  
//  }
//}

void phantomSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker slices: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
