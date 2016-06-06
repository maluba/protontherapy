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
/*
Author: Susanna Guatelli
*/
//
//    **********************************
//    *                                *
//    *     MedPhysicsList.cc          *
//    *                                *
//    **********************************
//
#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4HadronElasticPhysics.hh"      // MCS
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
//#include "G4UHadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonElasticPhysics.hh"

#include "G4IonBinaryCascadePhysics.hh"    // ion interactions
#include "G4IonPhysics.hh"
#include "G4IonQMDPhysics.hh"
#include "G4PrecoProtonBuilder.hh"

#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4AutoDelete.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "MedPhysicsList.hh"		// MedPhys header file
#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SpecialCuts.hh"

MedPhysicsList::MedPhysicsList():  G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
    defaultCutValue = 5.*mm;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;

  SetVerboseLevel(1); 
 
  // EM physics: 3 alternatives

  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Alternatively you can substitute this physics list
  // with the LowEnergy Livermore or LowEnergy Penelope: 
  // emPhysicsList = new G4EmLivermorePhysics();
  // Low Energy based on Livermore Evaluated Data Libraries
  //
  // Penelope physics
  //emPhysicsList = new G4EmPenelopePhysics();

  // Hadronic physics
  //
  hadronPhys.push_back( new G4DecayPhysics());
  hadronPhys.push_back( new G4RadioactiveDecayPhysics());
  hadronPhys.push_back( new G4IonBinaryCascadePhysics());
//  hadronPhys.push_back( new G4PrecoProtonBuilder());
  hadronPhys.push_back( new G4EmExtraPhysics());
  hadronPhys.push_back( new G4HadronElasticPhysics());
  hadronPhys.push_back( new G4HadronHElasticPhysics());
//  hadronPhys.push_back( new G4UHadronElasticPhysics());   //added
  hadronPhys.push_back( new G4StoppingPhysics());
  hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());

  // Decay physics
  decPhysicsList = new G4DecayPhysics();
  radDecayPhysicsList = new G4RadioactiveDecayPhysics();

}

MedPhysicsList::~MedPhysicsList()
{  
  delete decPhysicsList;
  delete radDecayPhysicsList;
  delete emPhysicsList;
  hadronPhys.clear();
  for(size_t i=0; i<hadronPhys.size(); i++)
  {
      delete hadronPhys[i];
  }
}

void MedPhysicsList::ConstructParticle()
{
    decPhysicsList -> ConstructParticle();
}

void MedPhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

  // EM physics
  //
  emPhysicsList -> ConstructProcess();
  em_config.AddModels();

  // Hadronic physics
  //
  for(size_t i=0; i < hadronPhys.size(); i++)
  {
        hadronPhys[i] -> ConstructProcess();
  }

  // decay physics list
   decPhysicsList -> ConstructProcess();
   radDecayPhysicsList -> ConstructProcess();

   // step limitation (as a full process)
   //
    AddStepMax();

// must I include the sensitive detector here?
}

/////////////////////////////////////////////////////////////////////////////
void MedPhysicsList::AddPhysicsList(const G4String& name)
{
    if (verboseLevel>1) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;
    
    ///////////////////////////////////
    //   ELECTROMAGNETIC MODELS
    ///////////////////////////////////

        if (name == "standard_opt3") {
        emName = name;
        delete emPhysicsList;
        hadronPhys.clear();
        emPhysicsList = new G4EmStandardPhysics_option3();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;

    ///////////////////////////////////
    //   HADRONIC MODELS
    ///////////////////////////////////

//    } else if (name == "local_ion_ion_inelastic") {
//        hadronPhys.push_back(new LocalIonIonInelasticPhysic());
//        locIonIonInelasticIsRegistered = true;
     }
     else if (name == "MedPhys_1") {
        
        // The MedPhys_1 physics list consists of:
        // StandardEM_opt3 for EM processes 
	// 20 bins/decade EM process initialization table
        // I = 75 eV excitation potential for water
        // G4UHadronElasticProcess, together with G4HadronElasticModel
	// PRECO, precompound model in place of G4IonBinaryCascadePhysics() for ion-ion inelastic 		interactions
        // --> The G4RadioactiveDecayPhysics is added
        
        AddPhysicsList("standard_opt3");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
     //	hadronPhys.push_back(new G4PrecoProtonBuilder());
     //   hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysics());
     //   hadronPhys.push_back( new G4UHadronElasticPhysics());
	hadronPhys.push_back( new G4HadronHElasticPhysics());
        hadronPhys.push_back( new G4StoppingPhysics());
     //   hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        
        G4cout << "MedPhys_1 PHYSICS LIST has been activated" << G4endl;
    
        
    } else if (name == "MedPhys_2") {
        
        // Same as MedPhys_1, except Lewis MCS model is replaced with Goudsmitsounderson
	// model of e- scattering
        
        AddPhysicsList("standard_opt3");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
     //   hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysics());
	hadronPhys.push_back( new G4HadronHElasticPhysics());
        hadronPhys.push_back( new G4StoppingPhysics());
     //   hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        
        G4cout << "MedPhys_2 PHYSICS LIST has been acivated" << G4endl;
        
    } 
    else if (name == "MedPhys_3") {
	//
	//BinaryCascadePhysics is used instead of PRECO for ion inelastic interaction
	//
      AddPhysicsList("QGSP_BIC_EMY");
      AddPhysicsList("standard_opt3");
      //emPhysicsList = new G4EmLivermorePhysics();
      hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
      hadronPhys.push_back( new G4EmExtraPhysics());
      hadronPhys.push_back( new G4HadronElasticPhysics());
     // hadronPhys.push_back( new G4UHadronElasticPhysics());       //added
      hadronPhys.push_back( new G4StoppingPhysics());
      hadronPhys.push_back( new G4IonBinaryCascadePhysics());
      hadronPhys.push_back( new G4NeutronTrackingCut());
      hadronPhys.push_back( new G4DecayPhysics());
    }   
     else {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
        << " is not defined"
        << G4endl;
    }
    
}

void MedPhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();

  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
        pmanager -> AddProcess(new G4StepLimiter(),  -1,-1,3); 
        //Numbers are process order-AlongStep AtRest PostStep (-1 means the process is not active)
  }

}


void MedPhysicsList::SetCuts()
{

if (verboseLevel >0){
        G4cout << "MedPhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

  // Definition of  threshold of production 
  // of secondary particles
  // This is defined in range.
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

// By default the low energy limit to produce 
// secondary particles is 990 eV.
// This value is correct when using the EM Standard Physics.
// When using the Low Energy Livermore this value can be 
// changed to 250 eV corresponding to the limit
// of validity of the physics models.
// Comment out following three lines if the 
// Standard electromagnetic Package is adopted.

//  G4double lowLimit = 250. * eV;
//  G4double highLimit = 100. * GeV;

//  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

  // Print the cuts 
  if (verboseLevel>0) DumpCutValuesTable();
}
