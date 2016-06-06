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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include "HSteppingAction.hh"
#include "HRunAction.hh"
#include "HAnalysis.hh"
#include "g4root.hh"

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
//#include "HistoManager.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

  using namespace std;
  ofstream myfile("myfilename");
// for the purpose of writing to a file , don't forget to put 
//#include <fstream> using namespace std; extern ofstream myfile; 
//in other source file where rest of data to be output is coming from

HSteppingAction::HSteppingAction(AnalysisManager* pAnalysis)
: G4UserSteppingAction()
{ 
analysis = pAnalysis;
fSecondary = 0;
}

HSteppingAction::~HSteppingAction()
{ 
}

//                       Sugi's suggestion                              
void HSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* aTrack = aStep -> GetTrack();
  
  //check if alive
  if(aTrack -> GetTrackStatus() == fAlive) {return ;}
  
  //Using get methods to access information
  
     G4int trackID            = aTrack->GetTrackID();
//     G4int stepNum = aStep->GetStepNumber(); // does the stepping object have this method?
     G4StepPoint* preStepPoint        = aStep->GetPreStepPoint();
     G4StepPoint* postStepPoint       = aStep->GetPostStepPoint();
     G4ParticleDefinition* particle   = aTrack->GetDefinition();
     G4ThreeVector coordOfPreStep     = aStep->GetPreStepPoint()->GetPosition();
     G4ThreeVector coordOfPostStep    = aStep->GetPostStepPoint()->GetPosition();
     G4double energyAtPreStep         = aStep->GetPreStepPoint()->GetKineticEnergy();
     G4double energyDeposit  = aStep->GetTotalEnergyDeposit();
     
     G4String creatorProcessName = "---";    // Get creation process      
     if (aTrack->GetCreatorProcess() != NULL){creatorProcessName = aTrack->GetCreatorProcess()->GetProcessName();}
     G4String createProcessVolume = aTrack->GetLogicalVolumeAtVertex()->GetName();
         
     G4TouchableHandle touchpreStep = preStepPoint->GetTouchableHandle();
     G4VPhysicalVolume* PreStepVolume = touchpreStep->GetVolume();
     G4String NamePreStepVolume = PreStepVolume->GetName();
     G4String particleName    = aTrack->GetDefinition()->GetParticleName();

 // Only record events inside the Phantom volume
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     if ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "RingPhys")  // replace your volume name here
         &&(aStep -> GetPostStepPoint() -> GetProcessDefinedStep() != NULL))
    {
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Primary Proton%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (particleName == "proton" && aTrack->GetParentID()  == 1)  // this is for primary proton energy deposit
              {
//#ifdef ANALYSIS_USE_ROOT  

          analysis -> FillPrimaryProtonHistogram(energyDeposit/keV); //store primary proton energy in a 1D histogram

//#endif
//              }

/*
//----------------------write to the --------------------------------------------------------
//{
  std::ofstream myfile ("primaryProtonEdep.txt", ios::out | ios::app);
  if (myfile.is_open())
  {
//   myfile << "Writing this to a file.\n"<< G4endl;
//   myfile << std::left << std::setw(20) << "Particle Name" << std::setw(20) << "x" << std::setw(20) << "y" 
//          << std::setw(20) << "z"<< std::setw(20) << "Kinetic Energy" 
//          << std::setw(20) << "Energy Deposited" << G4endl;

  G4double x = coordOfPostStep.x()/cm;
  G4double y = coordOfPostStep.y()/cm;
  G4double z = coordOfPostStep.z()/cm;

  G4double r = std::sqrt(x*x + y*y);

   myfile << std::left << std::setw(20) << particleName << std::setw(20) << x << std::setw(20) << y << std::setw(20) << z 
          << std::setw(20) << energyAtPreStep/MeV << std::setw(20) << r << std::setw(20) << energyDeposit/MeV << G4endl;

   myfile.close();
  }
*/
} 

//        {       
//         G4cout << std::left  << std::setw(9) << trackID << /*std::setw(9) << stepNum <<*/ std::setw(15) << particleName  << std::setw(19)
//              << creatorProcessName <</* std::setw(19) << postProcessName<<*/ std::setw(15) << std::setprecision(7) << energyAtPreStep/MeV
//          << std::setw(15) << coordOfPostStep.z()/cm<< std::setw(15) <<energyDeposit/MeV  <<G4endl; 
//               // check your beam direction and replace Z()
//        }


       // Retrieve information about the secondaries originated in the phantom
//          G4SteppingManager*  steppingManager = fpSteppingManager;
          G4TrackVector* fSecondary = steppingManager -> GetfSecondary();

          for(size_t lp1 = 0; lp1 < (*fSecondary).size(); lp1++)
          { 
             G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); //volume where secondary originates
             G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  
             G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy();
             G4ThreeVector direction = (*fSecondary)[lp1]->GetMomentumDirection(); 
             G4ThreeVector coordinate  = (*fSecondary)[lp1]->GetPosition();   
//             G4double SecondaryenergyDeposit = (*fSecondary)[lp1]->GetTotalEnergyDeposit();  
             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Secondary Protons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
          if (secondaryParticleName == "proton" ) // this is for secondary proton energy deposit
           {
#ifdef ANALYSIS_USE_ROOT
//          G4double energy = (*fSecondary)[lp1]  -> GetKineticEnergy();
//          G4double energyDeposit  = aStep->GetTotalEnergyDeposit();
//           Store the initial energy of secondary particles in a 1D histogram
          analysis -> FillScondaryProtonHistogram(energyDeposit/keV);
#endif
//           }
/*
//{
// writing to file
  std::ofstream myfile ("secondaryProtonEdep.txt", ios::out | ios::app);
  if (myfile.is_open())
  {
//   myfile << "Writing this to a file.\n"<< G4endl;
//   myfile << std::left << std::setw(20) << "Particle Name" << std::setw(20) << "x" << std::setw(20) << "y" 
//          << std::setw(20) << "z"<< std::setw(20) << "Kinetic Energy" 
//          << std::setw(20) << "Energy Deposited" << G4endl;

   myfile << std::left << std::setw(20) << secondaryParticleName << std::setw(20) << coordinate.x()/cm << std::setw(20) 
          << coordinate.y()/cm << std::setw(20) << coordinate.z()/cm 
          << std::setw(20) << secondaryParticleKineticEnergy/MeV << std::setw(20) << energyDeposit/MeV << G4endl;

   myfile.close();
  }
*/
} 

//              {       
//           G4cout << std::left  << std::setw(15) << particleName  << std::setw(19)  
//                       << coordOfPostStep.z()/cm<< std::setw(15) <<energyDeposit/MeV  <<G4endl;
//                      // check your beam direction and replace Z()
//               }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Recoil ions & other Secondaries%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          else if((secondaryParticleName == "neutron")||
                 (secondaryParticleName == "alpha") ||
                 (secondaryParticleName == "deuteron") || 
                 (secondaryParticleName == "triton") || 
                 (secondaryParticleName == "He3") || 
//               (secondaryParticleName == "O16") ||
	         (secondaryParticleName =="GenericIon"))
           {
#ifdef ANALYSIS_USE_ROOT
//          G4double energy = (*fSecondary)[lp1]  -> GetKineticEnergy();
//          G4double energyDeposit  = aStep->GetTotalEnergyDeposit();
//           Store the initial energy of secondary particles in a 1D histogram
          analysis -> FillOtherSecondariesHistogram(energyDeposit/keV);
#endif
//           }
/*
//{ //shutdonw the printing to a file
       // writing to file
       std::ofstream myfile ("otherSecondaryParticleEdep.txt", ios::out | ios::app | ios::binary);
       if (myfile.is_open())
       {
//   myfile << "Writing this to a file.\n"<< G4endl;
//   myfile << std::left << std::setw(20) << "Particle Name" << std::setw(20) << "x" << std::setw(20) << "y" 
//          << std::setw(20) << "z"<< std::setw(20) << "Kinetic Energy" 
//          << std::setw(20) << "Energy Deposited" << G4endl;

       myfile << std::left << std::setw(20) << secondaryParticleName << std::setw(20) << coordinate.x()/cm << std::setw(20) 
              << coordinate.y()/cm << std::setw(20) << coordinate.z()/cm 
              << std::setw(20) << secondaryParticleKineticEnergy/MeV << std::setw(20) << energyDeposit/MeV << G4endl;

      myfile.close();
      }
  */
} 
          }

    }

}


