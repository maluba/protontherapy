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
// $Id: exampleB2a.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file exampleB2a.cc
/// \brief Main program of the B2a example

// This program uses the SensitiveDetector, hits and user hooks (G4UserRunAction and G4UserEventAction for run summary)

#include "HDetectorConstruction.hh"
#include "HActionInitialization.hh"
//#include "AnalysisManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsListHelper.hh"
#include "G4VUserPhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "MedPhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
//#include "QGSP_BERT_EMY.hh" //"FTFP_BERT.hh" This defines the physicsList in this code
#include "G4StepLimiterPhysics.hh"

//#include "HRunAction.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " hadrontherapy [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 6;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  


  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  G4RunManager * runManager = new G4RunManager;
#endif

//  AnalysisManager* analysis = new AnalysisManager();

  // Set mandatory initialization classes
  //
  HDetectorConstruction* detConstruction = new HDetectorConstruction();
  runManager->SetUserInitialization(detConstruction);
  
//initialize the physicsList
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  G4String physname = "";
  
  //method of calling a physicslist
  char* path = getenv("PHYSLIST");
  if(path){physname = G4String(path); }
  
  if(physname != "" && factory.IsReferencePhysList(physname))
  {
    phys = factory.GetReferencePhysList(physname);
  }
//  if(phys)
//  {
//    G4cout <<"Going to register the phantom region" << G4endl;
//    phys->RegisterPhysics(new )
//  }
//  else
//  {
//    G4cout <<"Use the user PhysicsList" << G4endl;
//    phys = new QGSP_BERT;
//  }
  else
  {
    G4cout << "Using MedPhysicsList()" << G4endl;
    //phys = new MedPhysicsList();
  }



  phys->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(phys);
  
//  G4VModularPhysicsList* physicsList = new QGSP_BERT; //QGSP_BIC_EMY;  //new QGSP_BERT;
//  runManager->SetUserInitialization(physicsList);
  
  //Set User action classes
  HActionInitialization* actionInitialization
     = new HActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !


#ifdef G4VIS_USE
  delete visManager;
//  delete analysis;
#endif

  delete runManager;
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
