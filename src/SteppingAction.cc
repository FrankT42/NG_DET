#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include <fstream>
#include "ActionInitialization.hh"
#include <sstream>

#include "G4Step.hh"

//======================================================
SteppingAction::SteppingAction(DetectorConstruction* det,
			       EventAction* evt)
  :G4UserSteppingAction(),
   fDetector(det), fEventAction(evt)
{
  Tout.open("./Outputs/Labels.dat");
}
//=====================================================
SteppingAction::~SteppingAction()
{}
//=====================================================
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  std::stringstream ss;
  //std::ofstream Tout("./Outputs/Types.dat");
  G4Track* track = aStep->GetTrack();
  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();
  if(ParticleName=="neutron"||ParticleName=="gamma"){
    ss<<ParticleName<<std::endl;
    Tout<<ss.str();
  }
  G4VPhysicalVolume* volume
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

 G4double GlobalTime = track->GetGlobalTime();
 G4double DeltaTime = aStep->GetDeltaTime();
 G4ThreeVector location = track->GetPosition();


//collect energy ect step by step
  if(volume ==fDetector->GetSiPM()){
    if(ParticleName == "opticalphoton"){
      fEventAction->AddSiPM(1);
      fEventAction->SiPMTrack(GlobalTime, DeltaTime, location);
      track->SetTrackStatus(fStopAndKill);
    }
    else{
      return;
    }
  }
  if( ParticleName=="opticalphoton"){return;}
 G4double edep = aStep->GetTotalEnergyDeposit();

 if(volume == fDetector->GetScint()) fEventAction->AddEdep(edep);
 


 const std::vector<const G4Track*>* secondaries =
                                            aStep->GetSecondaryInCurrentStep();

 if (secondaries->size()>0) {
   for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition()){
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
		  == "Scintillation")fEventAction->AddScint(1);
 
	   }
	}
   
 }
 }
 
 }
