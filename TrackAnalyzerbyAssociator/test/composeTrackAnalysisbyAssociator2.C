#include <string>
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"

#define debug 1

#ifndef PI
#define PI 3.1415926535
#endif

#define PtScale 20
#define PtNumber 11

void composeTrackAnalysisbyAssociator2(string ParaListName, int ParaNumber) {

    if(debug) cout << ParaListName << endl;
    string PtListName;
    ifstream composeParaList;
    composeParaList.open(ParaListName.c_str());

    string OutputPlotNamepreFix = ParaListName + "_";
    string OutputPlotNameFix = ".png";
    string SaveName;

    unsigned int EventNumber;
    unsigned int trackingParticleMatch;
    double recTrackPurity;
    double recTrackrefMomentum;
    double recTrackrefPhi;
    double recTrackrefEta;
    double recTrackinnerMomentum;
    double recTrackinnerPhi;
    double recTrackinnerEta;
    unsigned int recTrackinnerValid;
    double recTrackouterMomentum;
    double recTrackouterPhi;
    double recTrackouterEta;
    unsigned int recTrackouterValid;
    double simTrackinnerMomentum;
    double simTrackinnerPhi;
    double simTrackinnerEta;
    unsigned int simTrackinnerMatch;
    double simTrackouterMomentum;
    double simTrackouterPhi;
    double simTrackouterEta;
    unsigned int simTrackouterMatch;
    double recTrackinnerMomentumofTSOS;
    double recTrackinnerPhiofTSOS;
    double recTrackinnerEtaofTSOS;
    unsigned int recTrackinnerValidofTSOS;
    double recTrackouterMomentumofTSOS;
    double recTrackouterPhiofTSOS;
    double recTrackouterEtaofTSOS;
    unsigned int recTrackouterValidofTSOS;
    double recTrackimpactMomentumofTSOS;
    double recTrackimpactPhiofTSOS;
    double recTrackimpactEtaofTSOS;
    unsigned int recTrackimpactValidofTSOS;
    int recTrackCharge;
    double simTrackMomentumPt;
    double simTrackPhi;
    double simTrackEta;
    int simTrackCharge;

    TObjArray* myEfficiencyHist = new TObjArray();
    TObjArray* myParticleHist = new TObjArray();
    TObjArray* mySTAHist = new TObjArray();
    TObjArray* myChargeCheckHist = new TObjArray();
    TObjArray* myDeltaPtHist = new TObjArray();
    TObjArray* myDeltaPhiHist = new TObjArray();
    TObjArray* myDeltaEtaHist = new TObjArray();

    vector<string> ParaName;
    ParaName.clear();

    vector<string> PtName;
    PtName.clear();

    for(int ParaIndex = 0; ParaIndex < ParaNumber; ParaIndex++) {
        getline(composeParaList, PtListName);
        ParaName.push_back(PtListName);

        string TempHistName;
        TempHistName = PtListName + "_Efficiency2simPt";
        TH1D* Efficiency2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_Particle2simPt";
        TH1D* Particle2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_STA2simPt";
        TH1D* STA2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_InverseChargeRato2simPt";
        TH1D* InverseChargeRato2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_DeltaPt2simPt";
        TH1D* DeltaPt2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_DeltaPhi2simPt";
        TH1D* DeltaPhi2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_DeltaEta2simPt";
        TH1D* DeltaEta2simPtHist = new TH1D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber);
        TempHistName = PtListName + "_MaxPurity2simPt";
        TH2D* MaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 6, 0., 1.2);
        TempHistName = PtListName + "_Multiplicity2simPt";
        TH2D* Multiplicity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 10, 0., 10.);
        TempHistName = PtListName + "_ChargeCheck2simPt";
        TH2D* ChargeCheck2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 5, -2.5, 2.5);
        TempHistName = PtListName + "_simTrackMomentumPtmaxPurity2simPt";
        TH2D* simTrackMomentumPtmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, (int)5*PtScale, 0, PtScale);
        TempHistName = PtListName + "_simTrackPhimaxPurity2simPt";
        TH2D* simTrackPhimaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 314, -PI, PI);
        TempHistName = PtListName + "_simTrackEtamaxPurity2simPt";
        TH2D* simTrackEtamaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 400, -2.0, 2.0);
        TempHistName = PtListName + "_recTrackimpactMomentumofTSOSmaxPurity2simPt";
        TH2D* recTrackimpactMomentumofTSOSmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, (int)5*PtScale, 0, PtScale);
        TempHistName = PtListName + "_recTrackimpactPhiofTSOSmaxPurity2simPt";
        TH2D* recTrackimpactPhiofTSOSmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 314, -PI, PI);
        TempHistName = PtListName + "_recTrackimpactEtaofTSOSmaxPurity2simPt";
        TH2D* recTrackimpactEtaofTSOSmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 600, -3.0, 3.0);
        TempHistName = PtListName + "_recTrackimpactValidofTSOSmaxPurity2simPt";
        TH2D* recTrackimpactValidofTSOSmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(), PtNumber, 0, PtNumber, 2, 0., 2.);
        TempHistName = PtListName + "_DeltaPtmaxPurity2simPt";
        TH2D* DeltaPtmaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(),  PtNumber, 0, PtNumber, (int)5*PtScale, -1.*PtScale, PtScale);
        TempHistName = PtListName + "_DeltaPhimaxPurity2simPt";
        TH2D* DeltaPhimaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(),  PtNumber, 0, PtNumber, 314, -PI, PI);
        TempHistName = PtListName + "_DeltaEtamaxPurity2simPt";
        TH2D* DeltaEtamaxPurity2simPtHist = new TH2D(TempHistName.c_str(), TempHistName.c_str(),  PtNumber, 0, PtNumber, 400, -2.0, 2.0);

        ifstream composePtList;
        composePtList.open(PtListName.c_str());
        string FileName;
        int ParaNumber = 11;
        for(int PtIndex = 0; PtIndex < PtNumber; PtIndex++) {
            getline(composePtList, FileName);

            PtName.push_back(FileName);

            string fullFileName = "data/"+ FileName + ".root";
            if(debug) cout << FileName << endl;
            TFile* RootFile = TFile::Open(fullFileName.c_str());

            TTree* T1 = (TTree*)RootFile->Get("ExTree");
            T1->SetBranchAddress("EventNumber", &EventNumber);
            T1->SetBranchAddress("trackingParticleMatch", &trackingParticleMatch);
            T1->SetBranchAddress("recTrackPurity", &recTrackPurity);
            T1->SetBranchAddress("recTrackrefMomentum", &recTrackrefMomentum);
            T1->SetBranchAddress("recTrackrefPhi", &recTrackrefPhi);
            T1->SetBranchAddress("recTrackrefEta", &recTrackrefEta);
            T1->SetBranchAddress("recTrackinnerMomentum", &recTrackinnerMomentum);
            T1->SetBranchAddress("recTrackinnerPhi", &recTrackinnerPhi);
            T1->SetBranchAddress("recTrackinnerEta", &recTrackinnerEta);
            T1->SetBranchAddress("recTrackinnerValid", &recTrackinnerValid);
            T1->SetBranchAddress("recTrackouterMomentum", &recTrackouterMomentum);
            T1->SetBranchAddress("recTrackouterPhi", &recTrackouterPhi);
            T1->SetBranchAddress("recTrackouterEta", &recTrackouterEta);
            T1->SetBranchAddress("recTrackouterValid", &recTrackouterValid);
            T1->SetBranchAddress("simTrackinnerMomentum", &simTrackinnerMomentum);
            T1->SetBranchAddress("simTrackinnerPhi", &simTrackinnerPhi);
            T1->SetBranchAddress("simTrackinnerEta", &simTrackinnerEta);
            T1->SetBranchAddress("simTrackinnerMatch", &simTrackinnerMatch);
            T1->SetBranchAddress("simTrackouterMomentum", &simTrackouterMomentum);
            T1->SetBranchAddress("simTrackouterPhi", &simTrackouterPhi);
            T1->SetBranchAddress("simTrackouterEta", &simTrackouterEta);
            T1->SetBranchAddress("simTrackouterMatch", &simTrackouterMatch);
            T1->SetBranchAddress("recTrackinnerMomentumofTSOS", &recTrackinnerMomentumofTSOS);
            T1->SetBranchAddress("recTrackinnerPhiofTSOS", &recTrackinnerPhiofTSOS);
            T1->SetBranchAddress("recTrackinnerEtaofTSOS", &recTrackinnerEtaofTSOS);
            T1->SetBranchAddress("recTrackinnerValidofTSOS", &recTrackinnerValidofTSOS);
            T1->SetBranchAddress("recTrackouterMomentumofTSOS", &recTrackouterMomentumofTSOS);
            T1->SetBranchAddress("recTrackouterPhiofTSOS", &recTrackouterPhiofTSOS);
            T1->SetBranchAddress("recTrackouterEtaofTSOS", &recTrackouterEtaofTSOS);
            T1->SetBranchAddress("recTrackouterValidofTSOS", &recTrackouterValidofTSOS);
            T1->SetBranchAddress("recTrackimpactMomentumofTSOS", &recTrackimpactMomentumofTSOS);
            T1->SetBranchAddress("recTrackimpactPhiofTSOS", &recTrackimpactPhiofTSOS);
            T1->SetBranchAddress("recTrackimpactEtaofTSOS", &recTrackimpactEtaofTSOS);
            T1->SetBranchAddress("recTrackimpactValidofTSOS", &recTrackimpactValidofTSOS);
            T1->SetBranchAddress("recTrackCharge", &recTrackCharge);
            T1->SetBranchAddress("simTrackMomentumPt", &simTrackMomentumPt);
            T1->SetBranchAddress("simTrackPhi", &simTrackPhi);
            T1->SetBranchAddress("simTrackEta", &simTrackEta);
            T1->SetBranchAddress("simTrackCharge", &simTrackCharge);

            unsigned int trackingParticleMatch_temp;
            unsigned int efficiency_temp;
            double recTrackPurity_temp;
            double recTrackrefMomentum_temp;
            double recTrackrefPhi_temp;
            double recTrackrefEta_temp;
            double recTrackinnerMomentum_temp;
            double recTrackinnerPhi_temp;
            double recTrackinnerEta_temp;
            unsigned int recTrackinnerValid_temp;
            double recTrackouterMomentum_temp;
            double recTrackouterPhi_temp;
            double recTrackouterEta_temp;
            unsigned int recTrackouterValid_temp;
            double simTrackinnerMomentum_temp;
            double simTrackinnerPhi_temp;
            double simTrackinnerEta_temp;
            unsigned int simTrackinnerMatch_temp;
            double simTrackouterMomentum_temp;
            double simTrackouterPhi_temp;
            double simTrackouterEta_temp;
            unsigned int simTrackouterMatch_temp;
            double recTrackinnerMomentumofTSOS_temp;
            double recTrackinnerPhiofTSOS_temp;
            double recTrackinnerEtaofTSOS_temp;
            unsigned int recTrackinnerValidofTSOS_temp;
            double recTrackouterMomentumofTSOS_temp;
            double recTrackouterPhiofTSOS_temp;
            double recTrackouterEtaofTSOS_temp;
            unsigned int recTrackouterValidofTSOS_temp;
            double recTrackimpactMomentumofTSOS_temp;
            double recTrackimpactPhiofTSOS_temp;
            double recTrackimpactEtaofTSOS_temp;
            unsigned int recTrackimpactValidofTSOS_temp;
            int recTrackCharge_temp;
            double simTrackMomentumPt_temp;
            double simTrackPhi_temp;
            double simTrackEta_temp;
            int simTrackCharge_temp;

            int Nentries = T1->GetEntries(); 
            for(int i = 0; i < Nentries; i++) { 
                T1->GetEntry(i);
                if(trackingParticleMatch == 0) {
                    MaxPurity2simPtHist->Fill(simTrackMomentumPt, 0);
                    Multiplicity2simPtHist->Fill(simTrackMomentumPt, 0);
                    //int tempParticleBinNumber = Particle2simPtHist->FindBin(simTrackMomentumPt);
                    int tempParticleBinNumber = PtIndex + 1;
                    double tempParticleBinValue = Particle2simPtHist->GetBinContent(tempParticleBinNumber);
                    tempParticleBinValue += 1.;
                    Particle2simPtHist->SetBinContent(tempParticleBinNumber, tempParticleBinValue);
                }
                else {
                    efficiency_temp = 1;
                    trackingParticleMatch_temp = trackingParticleMatch;
                    recTrackPurity_temp = recTrackPurity;
                    recTrackrefMomentum_temp = recTrackrefMomentum;
                    recTrackrefPhi_temp = recTrackrefPhi;
                    recTrackrefEta_temp = recTrackrefEta;
                    recTrackinnerMomentum_temp = recTrackinnerMomentum;
                    recTrackinnerPhi_temp = recTrackinnerPhi;
                    recTrackinnerEta_temp = recTrackinnerEta;
                    recTrackinnerValid_temp = recTrackinnerValid;
                    recTrackouterMomentum_temp = recTrackouterMomentum;
                    recTrackouterPhi_temp = recTrackouterPhi;
                    recTrackouterEta_temp = recTrackouterEta;
                    recTrackouterValid_temp = recTrackouterValid;
                    simTrackinnerMomentum_temp = simTrackinnerMomentum;
                    simTrackinnerPhi_temp = simTrackinnerPhi;
                    simTrackinnerEta_temp = simTrackinnerEta;
                    simTrackinnerMatch_temp = simTrackinnerMatch;
                    simTrackouterMomentum_temp = simTrackouterMomentum;
                    simTrackouterPhi_temp = simTrackouterPhi;
                    simTrackouterEta_temp = simTrackouterEta;
                    simTrackouterMatch_temp = simTrackouterMatch;
                    recTrackinnerMomentumofTSOS_temp = recTrackinnerMomentumofTSOS;
                    recTrackinnerPhiofTSOS_temp = recTrackinnerPhiofTSOS;
                    recTrackinnerEtaofTSOS_temp = recTrackinnerEtaofTSOS;
                    recTrackinnerValidofTSOS_temp = recTrackinnerValidofTSOS;
                    recTrackouterMomentumofTSOS_temp = recTrackouterMomentumofTSOS;
                    recTrackouterPhiofTSOS_temp = recTrackouterPhiofTSOS;
                    recTrackouterEtaofTSOS_temp = recTrackouterEtaofTSOS;
                    recTrackouterValidofTSOS_temp = recTrackouterValidofTSOS;
                    recTrackimpactMomentumofTSOS_temp = recTrackimpactMomentumofTSOS;
                    recTrackimpactPhiofTSOS_temp = recTrackimpactPhiofTSOS;
                    recTrackimpactEtaofTSOS_temp = recTrackimpactEtaofTSOS;
                    recTrackimpactValidofTSOS_temp = recTrackimpactValidofTSOS;
                    recTrackCharge_temp = recTrackCharge;
                    simTrackMomentumPt_temp = simTrackMomentumPt;
                    simTrackPhi_temp = simTrackPhi;
                    simTrackEta_temp = simTrackEta;
                    simTrackCharge_temp = simTrackCharge;

                    bool nextStep = true;
                    while(nextStep) {
                        i++;
                        T1->GetEntry(i);
                        if(trackingParticleMatch <= trackingParticleMatch_temp)
                            nextStep = false;
                        else
                            trackingParticleMatch_temp = trackingParticleMatch;
                        if(nextStep == true && recTrackPurity_temp < recTrackPurity) {
                            if(debug) cout << "step another match, trackingParticleMatch_temp: " << trackingParticleMatch_temp << endl;
                            //trackingParticleMatch_temp = trackingParticleMatch;
                            recTrackPurity_temp = recTrackPurity;
                            recTrackrefMomentum_temp = recTrackrefMomentum;
                            recTrackrefPhi_temp = recTrackrefPhi;
                            recTrackrefEta_temp = recTrackrefEta;
                            recTrackinnerMomentum_temp = recTrackinnerMomentum;
                            recTrackinnerPhi_temp = recTrackinnerPhi;
                            recTrackinnerEta_temp = recTrackinnerEta;
                            recTrackinnerValid_temp = recTrackinnerValid;
                            recTrackouterMomentum_temp = recTrackouterMomentum;
                            recTrackouterPhi_temp = recTrackouterPhi;
                            recTrackouterEta_temp = recTrackouterEta;
                            recTrackouterValid_temp = recTrackouterValid;
                            simTrackinnerMomentum_temp = simTrackinnerMomentum;
                            simTrackinnerPhi_temp = simTrackinnerPhi;
                            simTrackinnerEta_temp = simTrackinnerEta;
                            simTrackinnerMatch_temp = simTrackinnerMatch;
                            simTrackouterMomentum_temp = simTrackouterMomentum;
                            simTrackouterPhi_temp = simTrackouterPhi;
                            simTrackouterEta_temp = simTrackouterEta;
                            simTrackouterMatch_temp = simTrackouterMatch;
                            recTrackinnerMomentumofTSOS_temp = recTrackinnerMomentumofTSOS;
                            recTrackinnerPhiofTSOS_temp = recTrackinnerPhiofTSOS;
                            recTrackinnerEtaofTSOS_temp = recTrackinnerEtaofTSOS;
                            recTrackinnerValidofTSOS_temp = recTrackinnerValidofTSOS;
                            recTrackouterMomentumofTSOS_temp = recTrackouterMomentumofTSOS;
                            recTrackouterPhiofTSOS_temp = recTrackouterPhiofTSOS;
                            recTrackouterEtaofTSOS_temp = recTrackouterEtaofTSOS;
                            recTrackouterValidofTSOS_temp = recTrackouterValidofTSOS;
                            recTrackimpactMomentumofTSOS_temp = recTrackimpactMomentumofTSOS;
                            recTrackimpactPhiofTSOS_temp = recTrackimpactPhiofTSOS;
                            recTrackimpactEtaofTSOS_temp = recTrackimpactEtaofTSOS;
                            recTrackimpactValidofTSOS_temp = recTrackimpactValidofTSOS;
                            recTrackCharge_temp = recTrackCharge;
                            simTrackMomentumPt_temp = simTrackMomentumPt;
                            simTrackPhi_temp = simTrackPhi;
                            simTrackEta_temp = simTrackEta;
                            simTrackCharge_temp = simTrackCharge;
                        }
                    }
                    i--;
                    //if(debug) cout << "Filling Multiplicity " << trackingParticleMatch_temp << endl;
                    MaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackPurity_temp);
                    Multiplicity2simPtHist->Fill(simTrackMomentumPt_temp, trackingParticleMatch_temp);
                    ChargeCheck2simPtHist->Fill(simTrackMomentumPt_temp, simTrackCharge_temp*recTrackCharge_temp);
                    simTrackMomentumPtmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, simTrackMomentumPt_temp);
                    simTrackPhimaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, simTrackPhi_temp);
                    simTrackEtamaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, simTrackEta_temp);
                    recTrackimpactMomentumofTSOSmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactMomentumofTSOS_temp);
                    recTrackimpactPhiofTSOSmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactPhiofTSOS_temp);
                    recTrackimpactEtaofTSOSmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactEtaofTSOS_temp);
                    recTrackimpactValidofTSOSmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactValidofTSOS_temp);
                    DeltaPtmaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, (recTrackimpactMomentumofTSOS_temp-simTrackMomentumPt_temp)/simTrackMomentumPt_temp);
                    DeltaPhimaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactPhiofTSOS_temp-simTrackPhi_temp);
                    DeltaEtamaxPurity2simPtHist->Fill(simTrackMomentumPt_temp, recTrackimpactEtaofTSOS_temp-simTrackEta_temp);


                    //int tempParticleBinNumber = STA2simPtHist->FindBin(simTrackMomentumPt_temp);
                    int tempParticleBinNumber = PtIndex+1;
                    double tempParticleBinValue = Particle2simPtHist->GetBinContent(tempParticleBinNumber);
                    tempParticleBinValue += 1.;
                    Particle2simPtHist->SetBinContent(tempParticleBinNumber, tempParticleBinValue);
                    double tempSTABinValue = STA2simPtHist->GetBinContent(tempParticleBinNumber);
                    tempSTABinValue += 1.;                        
                    STA2simPtHist->SetBinContent(tempParticleBinNumber, tempSTABinValue);
                }
            }
            delete RootFile;
            delete T1;
        }

        for(int PtIndex = 1; PtIndex <= PtNumber; PtIndex++) {
            double ParticleBinValue = Particle2simPtHist->GetBinContent(PtIndex);
            double STABinValue = STA2simPtHist->GetBinContent(PtIndex);
            if(ParticleBinValue == 0.)
                ParticleBinValue += 1.;
            double EfficiencyBinValue = STABinValue / ParticleBinValue * 100.;
            double EfficiencyBinError = sqrt(EfficiencyBinValue * (100. - EfficiencyBinValue) / ParticleBinValue);
            cout << ParticleBinValue << ", " << STABinValue << ", " << EfficiencyBinValue << endl;
            Efficiency2simPtHist->SetBinContent(PtIndex, EfficiencyBinValue);
            Efficiency2simPtHist->SetBinError(PtIndex, EfficiencyBinError);

            TH1D* ChargeCheckHist = ChargeCheck2simPtHist->ProjectionY("ChargeCheck", PtIndex, PtIndex, "o");
            double ReverseChargeBinValue = ChargeCheckHist->GetBinContent(2);
            double CoverseChargeBinValue = ChargeCheckHist->GetBinContent(4);
            double TotalChargeBinValue = ReverseChargeBinValue + CoverseChargeBinValue;
            if(TotalChargeBinValue == 0.)
                TotalChargeBinValue += 1.;
            double ReverseChargeRato = ReverseChargeBinValue / TotalChargeBinValue;
            InverseChargeRato2simPtHist->SetBinContent(PtIndex, ReverseChargeRato);

            TH1D* DeltaPtHist = DeltaPtmaxPurity2simPtHist->ProjectionY("DeltaPt", PtIndex, PtIndex, "o");
            double DeltaPtMean = DeltaPtHist->GetMean();
            double DeltaPtRMS = DeltaPtHist->GetRMS();
            DeltaPt2simPtHist->SetBinContent(PtIndex, DeltaPtMean);
            DeltaPt2simPtHist->SetBinError(PtIndex, DeltaPtRMS);

            TH1D* DeltaPhiHist = DeltaPhimaxPurity2simPtHist->ProjectionY("DeltaPhi", PtIndex, PtIndex, "o");
            double DeltaPhiMean = DeltaPhiHist->GetMean();
            double DeltaPhiRMS = DeltaPhiHist->GetRMS();
            DeltaPhi2simPtHist->SetBinContent(PtIndex, DeltaPhiMean);
            DeltaPhi2simPtHist->SetBinError(PtIndex, DeltaPhiRMS);

            TH1D* DeltaEtaHist = DeltaEtamaxPurity2simPtHist->ProjectionY("DeltaEta", PtIndex, PtIndex, "o");
            double DeltaEtaMean = DeltaEtaHist->GetMean();
            double DeltaEtaRMS = DeltaEtaHist->GetRMS();
            DeltaEta2simPtHist->SetBinContent(PtIndex, DeltaEtaMean);
            DeltaEta2simPtHist->SetBinError(PtIndex, DeltaEtaRMS);
        }
        myEfficiencyHist->AddLast(Efficiency2simPtHist);
        myParticleHist->AddLast(Particle2simPtHist);
        mySTAHist->AddLast(STA2simPtHist);
        myChargeCheckHist->AddLast(InverseChargeRato2simPtHist);
        myDeltaPtHist->AddLast(DeltaPt2simPtHist);
    }

    double minX = 0;
    double minY = 0;
    double maxX = 110;
    double maxY = 40;

    TCanvas* myCanvas = new TCanvas("Canvas", "Canvas", 800, 600);
    myCanvas->cd();
    TPad* myPad = new TPad("Pad", "Pad", 0, 0, 1, 1);
    myPad->Draw();
    myPad->cd();

    ((TH1D*)(myParticleHist->At(0)))->SetStats(0);
    ((TH1D*)(myParticleHist->At(0)))->GetXaxis()->SetTitle("simPt/Gev");
    ((TH1D*)(myParticleHist->At(0)))->GetXaxis()->CenterTitle(1);
    ((TH1D*)(myParticleHist->At(0)))->Draw();
    for(int ParaIndex = 0; ParaIndex < ParaNumber; ParaIndex++) {
        ((TH1D*)(mySTAHist->At(ParaIndex)))->SetStats(0);
        ((TH1D*)(mySTAHist->At(ParaIndex)))->SetLineColor(kRed+ParaIndex);
        ((TH1D*)(mySTAHist->At(ParaIndex)))->Draw("same");
    }
    TLegend *STALeg = new TLegend(0.6,0.1,0.9,0.3);
    STALeg->SetBorderSize(1);
    TString LegKey = "ParticleTrack";
    STALeg->AddEntry(myParticleHist->At(0), LegKey, "lpf");
    for(int ParaIndex = 0; ParaIndex < ParaNumber; ParaIndex++) {
        LegKey = ParaName[ParaIndex];
        STALeg->AddEntry(mySTAHist->At(ParaIndex), LegKey, "lpf");
    }
    STALeg->Draw();
    SaveName = OutputPlotNamepreFix + "_STA2simPt" + OutputPlotNameFix;
    myCanvas->SaveAs(SaveName.c_str());

    myPad->Clear();
    myPad->Update();
    double YScale = myPad->GetUymax() / 110.;
    ((TH1D*)(myEfficiencyHist->At(0)))->GetXaxis()->SetTitle("simPt/Gev");
    ((TH1D*)(myEfficiencyHist->At(0)))->GetXaxis()->CenterTitle(1);
    ((TH1D*)(myEfficiencyHist->At(0)))->SetStats(0);
    ((TH1D*)(myEfficiencyHist->At(0)))->Scale(YScale);
    ((TH1D*)(myEfficiencyHist->At(0)))->SetLineColor(kRed);
    ((TH1D*)(myEfficiencyHist->At(0)))->Draw("same,ah");
    for(int ParaIndex = 1; ParaIndex < ParaNumber; ParaIndex++) {
        ((TH1D*)(myEfficiencyHist->At(ParaIndex)))->SetStats(0);
        ((TH1D*)(myEfficiencyHist->At(ParaIndex)))->Scale(YScale);
        ((TH1D*)(myEfficiencyHist->At(ParaIndex)))->SetLineColor(kRed+ParaIndex);
        ((TH1D*)(myEfficiencyHist->At(ParaIndex)))->Draw("same,ah");
    }
    myPad->Update();
    if(debug) cout << "Y: " << myPad->GetUymax() << endl;
    double YAxisMinValue=((TH1D*)(myEfficiencyHist->At(0)))->GetYaxis()->GetXmin();
    double YAxisMaxValue=((TH1D*)(myEfficiencyHist->At(0)))->GetYaxis()->GetXmax();
    int YAxisNBins=((TH1D*)(myEfficiencyHist->At(0)))->GetYaxis()->GetNbins();
    TGaxis* YAxis = new TGaxis(myPad->GetUxmin(), myPad->GetUymin(), myPad->GetUxmin(), myPad->GetUymax(), 0, 110, 510, "-R");
    YAxis->SetLineColor(kGreen);
    YAxis->SetLabelColor(kGreen);
    YAxis->SetTitle("Efficiency of STA for simPts");
    YAxis->CenterTitle(1);
    YAxis->Draw();
    double XAxisMinValue=((TH1D*)(myEfficiencyHist->At(0)))->GetXaxis()->GetXmin();
    double XAxisMaxValue=((TH1D*)(myEfficiencyHist->At(0)))->GetXaxis()->GetXmax();
    int XAxisNBins=((TH1D*)(myEfficiencyHist->At(0)))->GetXaxis()->GetNbins();
    TGaxis* XAxis = new TGaxis(myPad->GetUxmin(), myPad->GetUymin(), myPad->GetUxmax(), myPad->GetUymin(), XAxisMinValue, XAxisMaxValue, 510, "+L");
    XAxis->SetTitle("simPt/Gev");
    XAxis->CenterTitle(1);
    XAxis->Draw();
    TLegend *EffLeg = new TLegend(0.1,0.9,0.4,1.0);
    EffLeg->SetBorderSize(1);
    for(int ParaIndex = 0; ParaIndex < ParaNumber; ParaIndex++) {
        TString LegKey = ParaName[ParaIndex];
        EffLeg->AddEntry(myEfficiencyHist->At(ParaIndex), LegKey, "lpf");
    }
    EffLeg->Draw();
    SaveName = OutputPlotNamepreFix + "_Eff2simPt" + OutputPlotNameFix;
    myCanvas->SaveAs(SaveName.c_str());

    ((TH1D*)(myDeltaPtHist->At(0)))->SetStats(0);
    ((TH1D*)(myDeltaPtHist->At(0)))->GetXaxis()->SetTitle("simPt/Gev");
    ((TH1D*)(myDeltaPtHist->At(0)))->GetXaxis()->CenterTitle(1);
    ((TH1D*)(myDeltaPtHist->At(0)))->GetYaxis()->SetTitle("deltPt/simPt");
    ((TH1D*)(myDeltaPtHist->At(0)))->GetYaxis()->CenterTitle(1);
    ((TH1D*)(myDeltaPtHist->At(0)))->SetLineColor(kRed);
    ((TH1D*)(myDeltaPtHist->At(0)))->Draw("");
    for(int ParaIndex = 1; ParaIndex < ParaNumber; ParaIndex++) {
        ((TH1D*)(myDeltaPtHist->At(ParaIndex)))->SetStats(0);
        //((TH1D*)(myDeltaPtHist->At(ParaIndex)))->GetXaxis()->SetTitle("simPt/Gev");
        //((TH1D*)(myDeltaPtHist->At(ParaIndex)))->GetXaxis()->CenterTitle(1);
        //((TH1D*)(myDeltaPtHist->At(ParaIndex)))->GetYaxis()->SetTitle("deltPt/simPt");
        //((TH1D*)(myDeltaPtHist->At(ParaIndex)))->GetYaxis()->CenterTitle(1);
        ((TH1D*)(myDeltaPtHist->At(ParaIndex)))->SetLineColor(kRed+ParaIndex);
        ((TH1D*)(myDeltaPtHist->At(ParaIndex)))->Draw("same");
        //SaveName = OutputPlotNamepreFix + ParaName[ParaIndex] + "DeltaPt" + OutputPlotNameFix;
        //myCanvas->SaveAs(SaveName.c_str());
    }
    TLegend *PtLeg = new TLegend(0.6,0.8,0.9,0.9);
    PtLeg->SetBorderSize(1);
    for(int ParaIndex = 0; ParaIndex < ParaNumber; ParaIndex++) {
        TString LegKey = ParaName[ParaIndex];
        PtLeg->AddEntry(myDeltaPtHist->At(ParaIndex), LegKey, "lpf");
    }
    PtLeg->Draw();
    SaveName = OutputPlotNamepreFix + "_DeltaPt" + OutputPlotNameFix;
    myCanvas->SaveAs(SaveName.c_str());
}
