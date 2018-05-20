#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "g4root.hh"
#include <fstream>
#include <iostream>
//#include "G4VAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingAction.hh"
#include <sstream>
//==================================
HistoManager::HistoManager()
  :fFactoryOn(false)
{
  std::ofstream fout("./Outputs/Pulses.dat");
}
//==================================
HistoManager::~HistoManager()
{}
//==================================

void HistoManager::Book()
{
  //make the outputs
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");

  G4bool fileOpen = analysisManager->OpenFile("N_G_Det");
  if(!fileOpen){
    G4cerr<<"\n---->Cannot open "<<analysisManager->GetFileName()<<G4endl;
    return;
  }
  //make histograms
  analysisManager->CreateH1("Edep", "Edep in Scintillator (MeV)",300,0,6*MeV);
  analysisManager->CreateH1("Scint","Number of Scintillation Photons",200,0,100000);
  analysisManager->CreateH1("SiPM","Number of Photons on SiPM",200,0,10000);
  analysisManager->CreateH1("x", "xposition", 512, -15*mm, 15*mm);
  analysisManager->CreateH1("y", "yposition", 512, -15*mm, 15*mm);
  analysisManager->CreateH1("time", "time1", 2048, 0, 5*us);
  analysisManager->CreateH1("time2","time2", 1030, 0, 1*ns);
  analysisManager->CreateH2("location", "loc", 1024, -15*mm,15*mm,1024,-15*mm,15*mm);
  //make ntuples
  analysisManager->CreateNtuple("Ntuple1", "Edep");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Scint");
  analysisManager->CreateNtupleDColumn("SiPM");
  analysisManager->FinishNtuple();
  //Time and loc
  analysisManager->CreateNtuple("Ntuple2","Time");
  analysisManager->CreateNtupleDColumn("Time1");
  analysisManager->CreateNtupleDColumn("Time2");
  analysisManager->FinishNtuple();
  analysisManager->CreateNtuple("Ntuple3","TDC");
  analysisManager->CreateNtupleDColumn("Time_TDC");
  analysisManager->FinishNtuple();
 
  
  fFactoryOn = true;
  G4cout<<"\n---->Output File is Open in "<<analysisManager->GetFileName()<<"."
	<<analysisManager->GetFileType()<<G4endl;

  //  std::ofstream fout;
  //fout.open("./Outputs/Pulses.dat");
						
}

//=========================================
void HistoManager::Save()
{
  if(!fFactoryOn) return;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  // fout.close();
  G4cout<<"\n---->Histograms and ntuples saved\n"<<G4endl;
}
void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(ih,xbin,weight);
 
}
void HistoManager::FillNtuple(G4double Edep, G4int Scint,G4int SiPM_num){
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if(Edep>0){
  analysisManager->FillNtupleDColumn(0,0,Edep);
  analysisManager->FillNtupleDColumn(0,1,Scint);
  analysisManager->FillNtupleDColumn(0,2,SiPM_num);
  analysisManager->AddNtupleRow(0);
  }
}
void HistoManager::FillTimeAndLoc(G4double Edep,std::vector<G4double> x, std::vector<G4double> y, std::vector<G4double>t1,std::vector<G4double> t2, std::vector<G4double> t1_1)
{
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  std::stringstream ss;
  for(int i = 0; i<x.size();i++){
    analysisManager->FillH1(3,x[i],1);
    analysisManager->FillH1(4,y[i],1);
    analysisManager->FillH2(0,y[i],x[i],1);

  }
    for(int i = 0; i<t1.size();i++){
    analysisManager->FillH1(5,t1[i],1);
    if(t1[i]<10*us){
    analysisManager->FillNtupleDColumn(1,0,t1[i]);
    }
    analysisManager->FillH1(6,t2[i],1);
	if(t2[i]<10*us){
    analysisManager->FillNtupleDColumn(1,1,t2[i]);
    ss<<t1[i]<<" ";
	}
	 analysisManager->AddNtupleRow(1);
    // SteppingAction::GetStream()<<std::endl;
    

  }
    // auto d1 = min_element(t1.begin(),t1.end());
    //	auto d2 = min_element(t1_1.begin(),t1_1.end());
    if(Edep>0){
    double d1 = 1<<20;
    double d2 = 1<<20;
    for(int i=0;i<t1.size();i++){
      if(t1[i]<d1&&t1[i]>0){d1=t1[i];}
    }
        for(int i=0;i<t1_1.size();i++){
	  if(t1_1[i]<d2&&t1_1[i]>0){d1=t1_1[i];}
    }
    G4double dif = abs(d1-d2);
     if(dif>0){
	analysisManager->FillNtupleDColumn(2,0,dif*ns);
		analysisManager->AddNtupleRow(2);
		}
    }
     //	std::cout<<dif<<std::endl;
	ss<<std::endl;
    SteppingAction::GetStream()<<ss.str();
    //        for(int i = 0; i<t2.size();i++){
    // analysisManager->FillH1(6,t2[i],1);
    //analysisManager->FillNtupleDColumn(1,1,t2[i]);

    //}
    //	analysisManager->AddNtupleRow(1);
}
