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
// $Id: B4cCalorHit.cc 69586 2013-05-08 14:20:11Z gcosmo $
//
/// \file B4cCalorHit.cc
/// \brief Implementation of the B4cCalorHit class

#include "phantomHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

//THIS IS NECESSARY FOR MT MODE
G4ThreadLocal G4Allocator<phantomHit>* phantomHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomHit::phantomHit()
 : G4VHit(),
   fTrack(-1),
   fTrackID(-1),
   fStripNo(-1),
   fEdep(0.),
   fPos(G4ThreeVector()),
   fTrackLength(0.),
   fMomentum(G4ThreeVector())
//   fParticle(G4ParticleDefinition* pdef)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomHit::~phantomHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

phantomHit::phantomHit(const phantomHit& right)
  : G4VHit()
{
  fTrackID      = right.fTrackID;
  fEdep         = right.fEdep;
  fTrackLength  = right.fTrackLength;
  fStripNo    = right.fStripNo;
  fPos          = right.fPos;
  fMomentum    = right.fMomentum;
  fParticle    = right.fParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const phantomHit& phantomHit::operator=(const phantomHit& right)
{
  fTrackID      = right.fTrackID;
  fEdep         = right.fEdep;
  fTrackLength  = right.fTrackLength;
  fStripNo    = right.fStripNo;
  fPos          = right.fPos;
  fMomentum    = right.fMomentum;
  fParticle    = right.fParticle;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int phantomHit::operator==(const phantomHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void phantomHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " StripNo: " << fStripNo
     << "\tEdep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     <<"Position: "
     << std::setw(7) << G4BestUnit( fPos,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
