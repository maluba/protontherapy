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
// $Id: B4DetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
//
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class


#include "HDetectorConstruction.hh"
#include "phantomSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal //---->check if you really need the magnetic field
G4GlobalMagFieldMessenger* HDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HDetectorConstruction::HDetectorConstruction()
 : G4VUserDetectorConstruction(),fNbofRings(0), fNofSlices(0),
   fRingLV(NULL), //fSliceLV(NULL),
   fCheckOverlaps(true)

{
  fNbofRings = 400;
  fNofSlices = 300;
  fRingLV = new G4LogicalVolume*[fNbofRings];
//  fSliceLV = new G4LogicalVolume*[fNofSlices]; //an array holding slices along the Z-axis
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HDetectorConstruction::~HDetectorConstruction()
{
  delete [] fRingLV;
//  delete [] fSliceLV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HDetectorConstruction::DefineMaterials()
{
  // material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_WATER");
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double worldX = 2. *m;
  G4double worldY = 2. *m;
  G4double worldZ = 2. *m;


  G4double SliceThicknessZ = 1.*mm;//1. *mm; // deltaz


  // Get materials
  G4Material* worldMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* phantomMaterial = G4Material::GetMaterial("G4_WATER");

  if ( ! worldMaterial || ! phantomMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("HDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  G4VSolid* world
    = new G4Box("World",           // its name
                 worldX/2, worldY/2, worldZ/2); // its size

  G4LogicalVolume* worldLV
    = new G4LogicalVolume(world,           // its solid
                          worldMaterial,  // its material
                          "World");         // its name

  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        worldLV,          // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps


  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
//  fTrackerCylinderLV[copyNo]->SetVisAttributes(simpleBoxVisAtt);
  G4VisAttributes* TrackerCylVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  TrackerCylVisAtt->SetVisibility(true);

  //oooOO0OOooo.....oooOO0OOooo..Evan's Geometry, Alex's finish..oooOO0OOooo.....oooOO0OOooo

  G4ThreeVector FirstCylposition = G4ThreeVector(0,0,0);
//  G4double Radii[]  {10.*mm,7.*mm,5.*mm,3.*mm,1.*mm,10.*mm}; //varying radii of the rings

  G4double rmin, rmax = 0;          //r_0;

  for(G4int copyNo=0; copyNo<fNbofRings; copyNo++)
  {
    rmin = rmax;
    //rmax = rmax + 0.5 * mm; //Radii[copyNo];
    //making the last ring much bigger than the rest-------------
    if(rmin<fNbofRings*0.5-0.5)
    {
        rmax = rmax + 0.5 * mm; //Radii[copyNo];
    }
    else
        rmax = rmax + 5.5 * mm;
    //-----------------------------------------------------------

    G4VSolid* Rings = new G4Tubs("Ring_solid", rmin,rmax, SliceThicknessZ/2,0.*deg,360.*deg);

    fRingLV[copyNo] = new G4LogicalVolume(Rings, phantomMaterial,"RingLV", 0,0,0);

//   fRingLV[copyNo]->SetVisAttributes(simpleBoxVisAtt);

    //reproducing cylinders along the z-axis

    for(G4int CopyNo=0; CopyNo<fNofSlices; CopyNo++)
    {
        G4double FirstZposition = 0;
        G4double Zposition = FirstZposition + CopyNo * SliceThicknessZ;
        new G4PVPlacement(0, G4ThreeVector(0, 0, Zposition), fRingLV[copyNo], "RingPhys", worldLV, false,CopyNo/*, fCheckOverlaps*/);
    }

  }

  //
  // print parameters
  //
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The phantom is divided into " << fNbofRings << "\tconcentric rings of width" << SliceThicknessZ
         << "\n and replicated along the beam-axis (\t" << fNofSlices* SliceThicknessZ <<"\t of the z-axis)"
         << "\n------------------------------------------------------------\n";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.1 *mm; //*hz;
  fStepLimit = new G4UserLimits(maxStep);
  for(G4int CopyNo=0; CopyNo<fNofSlices; CopyNo++){
  fRingLV[CopyNo]->SetUserLimits(fStepLimit);
}
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //---------------------------------Rings------------------------------------
  G4String ringSDname = "/hadrontherapy/phantomSD";
  phantomSD* SDname1 = new phantomSD(ringSDname, "PhantomSliceHitsCollection", fNofSlices /*fNbofRings*/);
  SetSensitiveDetector("RingLV",SDname1,true);

  //one lv cannot have more than detector object, what if we swap them

  //
  // Magnetic field

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.

  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo
void HDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
