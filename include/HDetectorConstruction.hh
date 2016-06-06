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
// $Id: B4cDetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4cDetectorConstruction.hh
/// \brief Definition of the B4cDetectorConstruction class

#ifndef HDetectorConstruction_h
#define HDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"
//#include "HRunAction.hh"
//#include "AnalysisManager.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class HDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    HDetectorConstruction();
    virtual ~HDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
     
     // set methods
     void SetMaxStep (G4double );
     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    G4LogicalVolume** fRingLV;
    G4UserLimits* fStepLimit;            // pointer to user step limits
//    G4LogicalVolume** fSliceLV;  //redundant
    
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4int   fNofSlices;     // number of repica rings _|_ z-axis
    G4int   fNbofRings;    //number of rings ripling out
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

