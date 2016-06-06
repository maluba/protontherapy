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
// $Id: B4cCalorHit.hh 69223 2013-04-23 12:36:10Z gcosmo $
//
/// \file B4cCalorHit.hh
/// \brief Definition of the B4cCalorHit class

#ifndef phantomHit_h
#define phantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"

class phantomHit : public G4VHit
{
  public:
    phantomHit();
    phantomHit(const phantomHit&);
    virtual ~phantomHit();

    // operators
    const phantomHit& operator=(const phantomHit&);
    G4int operator==(const phantomHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    //void Add(G4double de, G4double dl);

    // Set methods
    void SetTrack  (G4int track)      { fTrack = track; };
    void SetTrackID  (G4int trackid)      { fTrackID = trackid; };
    void SetStripNo(G4int strip)      { fStripNo = strip; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
    void SetTrackLength (G4double tracklength) { fTrackLength = tracklength; };
    inline void SetMomentum(G4ThreeVector mom) { fMomentum = mom; };
    inline void SetParticle(G4ParticleDefinition* pdef) { fParticle = pdef; };

    inline void Add(G4double de, G4double dl) { 
      fEdep += de;
      fTrackLength += dl;
      }

    // Get methods
    inline G4int GetTrack() const     { return fTrack; };
    inline G4int GetTrackID() const     { return fTrackID; };
    inline G4int GetStripNo() const   { return fStripNo; };
    inline G4double GetEdep() const     { return fEdep; };
    inline G4ThreeVector GetPos() const { return fPos; };
    inline G4double GetTrackLength() const {return fTrackLength; };
    inline G4ThreeVector GetMomentum() { return fMomentum; };
    inline G4ParticleDefinition* GetParticle() { return fParticle; };

  private:
      G4int         fTrack;
      G4int         fTrackID;
      G4int         fStripNo;
      G4double      fEdep;
      G4ThreeVector fPos;
      G4ThreeVector fMomentum;
      G4double fTrackLength;
      G4ParticleDefinition* fParticle;
      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<phantomHit> phantomHitsCollection;

extern G4ThreadLocal G4Allocator<phantomHit>* phantomHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* phantomHit::operator new(size_t)
{
  if(!phantomHitAllocator)
    phantomHitAllocator = new G4Allocator<phantomHit>;
  void *hit;
  hit = (void *) phantomHitAllocator->MallocSingle();
  return hit;
}

inline void phantomHit::operator delete(void *hit)
{
  if(!phantomHitAllocator)
    phantomHitAllocator = new G4Allocator<phantomHit>;
  phantomHitAllocator->FreeSingle((phantomHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
