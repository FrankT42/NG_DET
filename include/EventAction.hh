#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
class RunAction;
class HistoManager;

class EventAction: public G4UserEventAction
{
public:
  EventAction(RunAction*, HistoManager*);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  inline void AddEdep(G4double de){fEdep+=de;};
  inline void AddScint(G4int dS){fScint++;};
  inline void AddSiPM(G4int dSi){fSiPM++;};
  void SiPMTrack(G4double global_time,G4double delta_time, G4ThreeVector location); 
private:
  RunAction* fRunAction;
  HistoManager* fHistoManager;

  G4double fEdep;
  G4int fSiPM;
  G4int fScint;
  G4int fPrintModulo;
  std::vector<G4double> x,y,t_g,t_d;

};
#endif
