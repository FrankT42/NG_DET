#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
//=============================================
//             Detector Simulation
//                 April 2018
//
//=============================================
//  Author:Frank Thomson
//  Simulation of Scintillation Detector
//  with neutron sensitivity and optical
//  transport
//=============================================

int main(int argc, char** argv){
  G4UIExecutive* ui = 0;

  if(argc == 1){
    ui = new G4UIExecutive(argc,argv);
  }
#ifdef G4MULTITHREADED
  G4int nThreads = 8;
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  //Lets get the mandatory stuff done
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization(det));

  //Make it SO!

  //Turn on the Telly (initialize visulization, but less posh)

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();//This shoudl be read in Brian Blessed's voice

  //Get the pointer to the UI, we will need that later
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if(!ui){//Damned thing is in batch mode, someone get a shovel
    G4String command = "/control/execute ";
    G4String fileName = argv[1]; //I know, my old way was cool, it will be back
    UImanager->ApplyCommand(command+fileName);
    
  }
  else{//Interactive mode. It must be too lonely in batch mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;//It was too good for this world :-(
  }
  //Termination...Find Sarah Connor
  delete visManager;
  delete runManager;

  return 0;//Skynet FTW

}
//=============================================
//           End of Main Method
//=============================================
