#include <string>
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TMath.h"

#define debug 0

#ifndef PI
#define PI 3.1415926535
#endif

#define PtScale 100

void TrackAnalysisbyAssociator(string FileName) {

    if(debug) cout << FileName << endl;
    TFile* RootFile = TFile::Open(FileName.c_str());

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

    TH1F* EfficiencyforTrackingParticle = new TH1F("EfficiencyforTrackingParticle", "EfficiencyforTrackingParticle", 2, 0, 2);
    TH1F* maxPurityperTrackingParticle = new TH1F("maxPurityperTrackingParticle", "maxPurityperTrackingParticle", 20, 0, 2);
    TH1F* MultiplicityperTrackingPartile = new TH1F("MultiplicityperTrackingPartile", "MultiplicityperTrackingPartile", 10, 0, 10);
    TH1F* recTrackrefMomentummaxPurity = new TH1F("recTrackrefMomentummaxPurity", "recTrackrefMomentummaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackrefPhimaxPurity = new TH1F("recTrackrefPhimaxPurity", "recTrackrefPhimaxPurity", 314, -PI, PI);
    TH1F* recTrackrefEtamaxPurity = new TH1F("recTrackrefEtamaxPurity", "recTrackrefEtamaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackinnerMomentummaxPurity = new TH1F("recTrackinnerMomentummaxPurity", "recTrackinnerMomentummaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackinnerPhimaxPurity = new TH1F("recTrackinnerPhimaxPurity", "recTrackinnerPhimaxPurity", 314, -PI, PI);
    TH1F* recTrackinnerEtamaxPurity = new TH1F("recTrackinnerEtamaxPurity", "recTrackinnerEtamaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackinnerValidmaxPurity = new TH1F("recTrackinnerValidmaxPurity", "recTrackinnerValidmaxPurity", 2, 0, 2);
    TH1F* recTrackouterMomentummaxPurity = new TH1F("recTrackouterMomentummaxPurity", "recTrackouterMomentummaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackouterPhimaxPurity = new TH1F("recTrackouterPhimaxPurity", "recTrackouterPhimaxPurity", 314, -PI, PI);
    TH1F* recTrackouterEtamaxPurity = new TH1F("recTrackouterEtamaxPurity", "recTrackouterEtamaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackouterValidmaxPurity = new TH1F("recTrackouterValidmaxPurity", "recTrackouterValidmaxPurity", 2, 0, 2);
    TH1F* simTrackMomentumPtmaxPurity = new TH1F("simTrackMomentumPtmaxPurity", "simTrackMomentumPtmaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* simTrackPhimaxPurity = new TH1F("simTrackPhimaxPurity", "simTrackPhimaxPurity", 314, -PI, PI);
    TH1F* simTrackEtamaxPurity = new TH1F("simTrackEtamaxPurity", "simTrackEtamaxPurity", 600, -3.0, 3.0);
    TH1F* simTrackinnerMomentummaxPurity = new TH1F("simTrackinnerMomentummaxPurity", "simTrackinnerMomentummaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* simTrackinnerPhimaxPurity = new TH1F("simTrackinnerPhimaxPurity", "simTrackinnerPhimaxPurity", 314, -PI, PI);
    TH1F* simTrackinnerEtamaxPurity = new TH1F("simTrackinnerEtamaxPurity", "simTrackinnerEtamaxPurity", 600, -3.0, 3.0);
    TH1F* simTrackinnerMatchmaxPurity = new TH1F("simTrackinnerMatchmaxPurity", "simTrackinnerMatchmaxPurity", 2, 0, 2);
    TH1F* simTrackouterMomentummaxPurity = new TH1F("simTrackouterMomentummaxPurity", "simTrackouterMomentummaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* simTrackouterPhimaxPurity = new TH1F("simTrackouterPhimaxPurity", "simTrackouterPhimaxPurity", 314, -PI, PI);
    TH1F* simTrackouterEtamaxPurity = new TH1F("simTrackouterEtamaxPurity", "simTrackouterEtamaxPurity", 600, -3.0, 3.0);
    TH1F* simTrackouterMatchmaxPurity = new TH1F("simTrackouterMatchmaxPurity", "simTrackouterMatchmaxPurity", 2, 0, 2);
    TH1F* recTrackinnerMomentumofTSOSmaxPurity = new TH1F("recTrackinnerMomentumofTSOSmaxPurity", "recTrackinnerMomentumofTSOSmaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackinnerPhiofTSOSmaxPurity = new TH1F("recTrackinnerPhiofTSOSmaxPurity", "recTrackinnerPhiofTSOSmaxPurity", 314, -PI, PI);
    TH1F* recTrackinnerEtaofTSOSmaxPurity = new TH1F("recTrackinnerEtaofTSOSmaxPurity", "recTrackinnerEtaofTSOSmaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackinnerValidofTSOSmaxPurity = new TH1F("recTrackinnerValidofTSOSmaxPurity", "recTrackinnerValidofTSOSmaxPurity", 2, 0, 2);
    TH1F* recTrackouterMomentumofTSOSmaxPurity = new TH1F("recTrackouterMomentumofTSOSmaxPurity", "recTrackouterMomentumofTSOSmaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackouterPhiofTSOSmaxPurity = new TH1F("recTrackouterPhiofTSOSmaxPurity", "recTrackouterPhiofTSOSmaxPurity", 314, -PI, PI);
    TH1F* recTrackouterEtaofTSOSmaxPurity = new TH1F("recTrackouterEtaofTSOSmaxPurity", "recTrackouterEtaofTSOSmaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackouterValidofTSOSmaxPurity = new TH1F("recTrackouterValidofTSOSmaxPurity", "recTrackouterValidofTSOSmaxPurity", 2, 0, 2);
    TH1F* recTrackimpactMomentumofTSOSmaxPurity = new TH1F("recTrackimpactMomentumofTSOSmaxPurity", "recTrackimpactMomentumofTSOSmaxPurity", (int)5*PtScale, 0, PtScale);
    TH1F* recTrackimpactPhiofTSOSmaxPurity = new TH1F("recTrackimpactPhiofTSOSmaxPurity", "recTrackimpactPhiofTSOSmaxPurity", 314, -PI, PI);
    TH1F* recTrackimpactEtaofTSOSmaxPurity = new TH1F("recTrackimpactEtaofTSOSmaxPurity", "recTrackimpactEtaofTSOSmaxPurity", 600, -3.0, 3.0);
    TH1F* recTrackimpactValidofTSOSmaxPurity = new TH1F("recTrackimpactValidofTSOSmaxPurity", "recTrackimpactValidofTSOSmaxPurity", 2, 0, 2);

    TH1F* deltaPtatimpactTSOSmaxPurity = new TH1F("deltaPtatimpactTSOSmaxPurity", "deltaPtatimpactTSOSmaxPurity", (int)10*PtScale, -1.*PtScale, PtScale);

    TH1D* Particle2simPtHist = new TH1D("Particle2simPt", "Particle2simPt", (int)(PtScale/2), 0, PtScale);
    TH1D* STA2simPtHist = new TH1D("STA2simPt", "STA2simPt", (int)(PtScale/2), 0, PtScale);

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
    double simTrackMomentumPt_temp;
    double simTrackPhi_temp;
    double simTrackEta_temp;

    int Nentries = T1->GetEntries(); 
    for(int i = 0; i < Nentries; i++) {
        T1->GetEntry(i);
        if(trackingParticleMatch == 0) {
            EfficiencyforTrackingParticle->Fill(0);
            maxPurityperTrackingParticle->Fill(0);
            MultiplicityperTrackingPartile->Fill(0);
            int tempParticleBinNumber = Particle2simPtHist->FindBin(simTrackMomentumPt);
            double tempParticleBinValue = Particle2simPtHist->GetBinContent(tempParticleBinNumber);
            tempParticleBinValue += 1.;
            Particle2simPtHist->SetBinContent(tempParticleBinNumber, tempParticleBinValue);
            //recTrackrefMomentummaxPurity->Fill(0);
            //recTrackrefPhimaxPurity->Fill(0);
            //recTrackrefEtamaxPurity->Fill(0);
            //recTrackinnerMomentummaxPurity->Fill(0);
            //recTrackinnerPhimaxPurity->Fill(0);
            //recTrackinnerEtamaxPurity->Fill(0);
            //recTrackinnerValidmaxPurity->Fill(0);
            //recTrackouterMomentummaxPurity->Fill(0);
            //recTrackouterPhimaxPurity->Fill(0);
            //recTrackouterEtamaxPurity->Fill(0);
            //recTrackouterValidmaxPurity->Fill(0);
            //simTrackMomentumPtmaxPurity->Fill(0);
            //simTrackPhimaxPurity->Fill(0);
            //simTrackEtamaxPurity->Fill(0);
            //simTrackinnerMomentummaxPurity->Fill(0);
            //simTrackinnerPhimaxPurity->Fill(0);
            //simTrackinnerEtamaxPurity->Fill(0);
            //simTrackinnerMatchmaxPurity->Fill(0);
            //simTrackouterMomentummaxPurity->Fill(0);
            //simTrackouterPhimaxPurity->Fill(0);
            //simTrackouterEtamaxPurity->Fill(0);
            //simTrackouterMatchmaxPurity->Fill(0);
            //recTrackinnerMomentumofTSOSmaxPurity->Fill(0);
            //recTrackinnerPhiofTSOSmaxPurity->Fill(0);
            //recTrackinnerEtaofTSOSmaxPurity->Fill(0);
            //recTrackinnerValidofTSOSmaxPurity->Fill(0);
            //recTrackouterMomentumofTSOSmaxPurity->Fill(0);
            //recTrackouterPhiofTSOSmaxPurity->Fill(0);
            //recTrackouterEtaofTSOSmaxPurity->Fill(0);
            //recTrackouterValidofTSOSmaxPurity->Fill(0);
            //recTrackimpactMomentumofTSOSmaxPurity->Fill(0);
            //recTrackimpactPhiofTSOSmaxPurity->Fill(0);
            //recTrackimpactEtaofTSOSmaxPurity->Fill(0);
            //recTrackimpactValidofTSOSmaxPurity->Fill(0);
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
            simTrackMomentumPt_temp = simTrackMomentumPt;
            simTrackPhi_temp = simTrackPhi;
            simTrackEta_temp = simTrackEta;

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
                    simTrackMomentumPt_temp = simTrackMomentumPt;
                    simTrackPhi_temp = simTrackPhi;
                    simTrackEta_temp = simTrackEta;
                }
            }
            i--;
            //if(debug) cout << "Filling Multiplicity " << trackingParticleMatch_temp << endl;
            EfficiencyforTrackingParticle->Fill(efficiency_temp);
            maxPurityperTrackingParticle->Fill(recTrackPurity_temp);
            MultiplicityperTrackingPartile->Fill(trackingParticleMatch_temp);
            recTrackrefMomentummaxPurity->Fill(recTrackrefMomentum);
            recTrackrefPhimaxPurity->Fill(recTrackrefPhi_temp);
            recTrackrefEtamaxPurity->Fill(recTrackrefEta_temp);
            recTrackinnerMomentummaxPurity->Fill(recTrackinnerMomentum_temp);
            recTrackinnerPhimaxPurity->Fill(recTrackinnerPhi_temp);
            recTrackinnerEtamaxPurity->Fill(recTrackinnerEta_temp);
            recTrackinnerValidmaxPurity->Fill(recTrackinnerValid_temp);
            recTrackouterMomentummaxPurity->Fill(recTrackouterMomentum_temp);
            recTrackouterPhimaxPurity->Fill(recTrackouterPhi_temp);
            recTrackouterEtamaxPurity->Fill(recTrackouterEta_temp);
            recTrackouterValidmaxPurity->Fill(recTrackouterValid_temp);
            simTrackMomentumPtmaxPurity->Fill(simTrackMomentumPt_temp);
            simTrackPhimaxPurity->Fill(simTrackPhi_temp);
            simTrackEtamaxPurity->Fill(simTrackEta_temp);
            simTrackinnerMomentummaxPurity->Fill(simTrackinnerMomentum_temp);
            simTrackinnerPhimaxPurity->Fill(simTrackinnerPhi_temp);
            simTrackinnerEtamaxPurity->Fill(simTrackinnerEta_temp);
            simTrackinnerMatchmaxPurity->Fill(simTrackinnerMatch_temp);
            simTrackouterMomentummaxPurity->Fill(simTrackouterMomentum_temp);
            simTrackouterPhimaxPurity->Fill(simTrackouterPhi_temp);
            simTrackouterEtamaxPurity->Fill(simTrackouterEta_temp);
            simTrackouterMatchmaxPurity->Fill(simTrackouterMatch_temp);
            recTrackinnerMomentumofTSOSmaxPurity->Fill(recTrackinnerMomentumofTSOS_temp);
            recTrackinnerPhiofTSOSmaxPurity->Fill(recTrackinnerPhiofTSOS_temp);
            recTrackinnerEtaofTSOSmaxPurity->Fill(recTrackinnerEtaofTSOS_temp);
            recTrackinnerValidofTSOSmaxPurity->Fill(recTrackinnerValidofTSOS_temp);
            recTrackouterMomentumofTSOSmaxPurity->Fill(recTrackouterMomentumofTSOS_temp);
            recTrackouterPhiofTSOSmaxPurity->Fill(recTrackouterPhiofTSOS_temp);
            recTrackouterEtaofTSOSmaxPurity->Fill(recTrackouterEtaofTSOS_temp);
            recTrackouterValidofTSOSmaxPurity->Fill(recTrackouterValidofTSOS_temp);
            recTrackimpactMomentumofTSOSmaxPurity->Fill(recTrackimpactMomentumofTSOS_temp);
            recTrackimpactPhiofTSOSmaxPurity->Fill(recTrackimpactPhiofTSOS_temp);
            recTrackimpactEtaofTSOSmaxPurity->Fill(recTrackimpactEtaofTSOS_temp);
            recTrackimpactValidofTSOSmaxPurity->Fill(recTrackimpactValidofTSOS_temp);
            
            deltaPtatimpactTSOSmaxPurity->Fill(recTrackimpactMomentumofTSOS_temp-simTrackMomentumPt_temp);
            int tempParticleBinNumber = STA2simPtHist->FindBin(simTrackMomentumPt_temp);
            double tempParticleBinValue = Particle2simPtHist->GetBinContent(tempParticleBinNumber);
            tempParticleBinValue += 1.;
            Particle2simPtHist->SetBinContent(tempParticleBinNumber, tempParticleBinValue);
            double tempSTABinValue = STA2simPtHist->GetBinContent(tempParticleBinNumber);
            tempSTABinValue += 1.;                        
            STA2simPtHist->SetBinContent(tempParticleBinNumber, tempSTABinValue);
        }
    }
	/*
    // get Eff2simPt
    for(int PtIndex = 1; PtIndex <= (int)(PtScale/2); PtIndex++) {
        double ParticleBinValue = Particle2simPt->GetBinContent(PtIndex);
        double STABinValue = STA2simPt->GetBinContent(PtIndex);
        if(ParticleBinValue == 0.)
            ParticleBinValue += 1.;
        double EfficiencyBinValue = STABinValue / ParticleBinValue * 100.;
        double EfficiencyBinError = sqrt(EfficiencyBinValue * (100. - EfficiencyBinValue) / ParticleBinValue);
        Eff2simPt->SetBinContent(PtIndex, EfficiencyBinValue);
        Eff2simPt->SetBinError(PtIndex, EfficiencyBinError);
    }
    */
    double minX = 0;
    double minY = 0;
    double maxX = 110;
    double maxY = 40;

    TCanvas* OutputCanvas = new TCanvas("Canvas", "Canvas", 800, 600);
    OutputCanvas->cd();
    TPad* myPad = new TPad("Pad", "Pad", 0, 0, 1, 1);
    myPad->Draw();
    myPad->cd();

    EfficiencyforTrackingParticle->Draw("");
    OutputCanvas->SaveAs("EfficiencyforTrackingParticle.png");
    myPad->SetLogy(1);
    maxPurityperTrackingParticle->Draw("");
    OutputCanvas->SaveAs("maxPurityperTrackingParticle.png");
    MultiplicityperTrackingPartile->Draw("");
    OutputCanvas->SaveAs("MultiplicityperTrackingPartile.png");
    myPad->SetLogy(0);
    recTrackrefMomentummaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackrefMomentummaxPurity.png");
    recTrackrefPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackrefPhimaxPurity.png");
    recTrackrefEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackrefEtamaxPurity.png");
    recTrackinnerMomentummaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerMomentummaxPurity.png");
    recTrackinnerPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerPhimaxPurity.png");
    recTrackinnerEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerEtamaxPurity.png");
    myPad->SetLogy(1);
    recTrackinnerValidmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerValidmaxPurity.png");
    myPad->SetLogy(0);
    recTrackouterMomentummaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterMomentummaxPurity.png");
    recTrackouterPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterPhimaxPurity.png");
    recTrackouterEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterEtamaxPurity.png");
    myPad->SetLogy(1);
    recTrackouterValidmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterValidmaxPurity.png");
    myPad->SetLogy(0);
    simTrackMomentumPtmaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackMomentumPtmaxPurity.png");
    simTrackPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackPhimaxPurity.png");
    simTrackEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackEtamaxPurity.png");
    simTrackinnerMomentummaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackinnerMomentummaxPurity.png");
    simTrackinnerPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackinnerPhimaxPurity.png");
    simTrackinnerEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackinnerEtamaxPurity.png");
    simTrackinnerMatchmaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackinnerMatchmaxPurity.png");
    simTrackouterMomentummaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackouterMomentummaxPurity.png");
    simTrackouterPhimaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackouterPhimaxPurity.png");
    simTrackouterEtamaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackouterEtamaxPurity.png");
    simTrackouterMatchmaxPurity->Draw("");
    OutputCanvas->SaveAs("simTrackouterMatchmaxPurity.png");
    recTrackinnerMomentumofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerMomentumofTSOSmaxPurity.png");
    recTrackinnerPhiofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerPhiofTSOSmaxPurity.png");
    recTrackinnerEtaofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerEtaofTSOSmaxPurity.png");
    myPad->SetLogy(1);
    recTrackinnerValidofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackinnerValidofTSOSmaxPurity.png");
    myPad->SetLogy(0);
    recTrackouterMomentumofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterMomentumofTSOSmaxPurity.png");
    recTrackouterPhiofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterPhiofTSOSmaxPurity.png");
    recTrackouterEtaofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterEtaofTSOSmaxPurity.png");
    myPad->SetLogy(1);
    recTrackouterValidofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackouterValidofTSOSmaxPurity.png");
    myPad->SetLogy(0);
    recTrackimpactMomentumofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackimpactMomentumofTSOSmaxPurity.png");
    recTrackimpactPhiofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackimpactPhiofTSOSmaxPurity.png");
    recTrackimpactEtaofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackimpactEtaofTSOSmaxPurity.png");
    recTrackimpactValidofTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("recTrackimpactValidofTSOSmaxPurity.png");
    deltaPtatimpactTSOSmaxPurity->Draw("");
    OutputCanvas->SaveAs("deltaPtatimpactTSOSmaxPurity.png");
}
