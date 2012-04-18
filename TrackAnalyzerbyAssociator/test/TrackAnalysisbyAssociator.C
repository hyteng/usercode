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
    
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat("");
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
    TH1F* recTrackEtamaxPurity = new TH1F("recTrackEtamaxPurity", "recTrackEtamaxPurity", 600, -3.0, 3.0);
    TH1F* EfficiencyEtamaxPurity = new TH1F("EfficiencyEtamaxPurity", "Efficiency w.r.t Eta", 600, -3.0, 3.0);

    TH1F* deltaPtatimpactTSOSmaxPurity = new TH1F("deltaPtatimpactTSOSmaxPurity", "deltaPtatimpactTSOSmaxPurity", (int)10*PtScale, -1.*PtScale, PtScale);

    TH1D* Particle2simPtHist = new TH1D("Particle2simPt", "Particle2simPt", (int)(PtScale/2), 0, PtScale);
    TH1D* STA2simPtHist = new TH1D("STA2simPt", "STA2simPt", (int)(PtScale/2), 0, PtScale);
    TH1D* Eff2simPtHist = new TH1D("STA2simPt", "STA2simPt", (int)(PtScale/2), 0, PtScale);
    TH1F* ChargeCheckHistTSOSmaxPurityHist = new TH1F("ChargeCheckHistTSOSmaxPurity", "ChargeCheckHistTSOSmaxPurity", 5, -2.5, 2.5);

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

        if(fabs(simTrackEta) >= 0.8)
            continue;

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
            simTrackMomentumPtmaxPurity->Fill(simTrackMomentumPt);
            simTrackPhimaxPurity->Fill(simTrackPhi);
            simTrackEtamaxPurity->Fill(simTrackEta);
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
                
                if(fabs(simTrackEta) >= 0.8)
                    continue;

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
            recTrackEtamaxPurity->Fill(simTrackEta_temp);
           
            deltaPtatimpactTSOSmaxPurity->Fill((recTrackimpactMomentumofTSOS_temp-simTrackMomentumPt_temp));
            int tempParticleBinNumber = STA2simPtHist->FindBin(simTrackMomentumPt_temp);
            double tempParticleBinValue = Particle2simPtHist->GetBinContent(tempParticleBinNumber);
            tempParticleBinValue += 1.;
            Particle2simPtHist->SetBinContent(tempParticleBinNumber, tempParticleBinValue);
            double tempSTABinValue = STA2simPtHist->GetBinContent(tempParticleBinNumber);
            tempSTABinValue += 1.;                        
            STA2simPtHist->SetBinContent(tempParticleBinNumber, tempSTABinValue);
            ChargeCheckHistTSOSmaxPurityHist->Fill(simTrackCharge*recTrackCharge);
        }
    }
	
    int BinNumber = simTrackEtamaxPurity->GetNbinsX();
    for(int i = 1; i <= BinNumber; i++) {
        double simTrackNumber = simTrackEtamaxPurity->GetBinContent(i);
        double recTrackNumber = recTrackEtamaxPurity->GetBinContent(i);
        if(simTrackNumber != 0.) {
            double BinEffMean = 100. * recTrackNumber / simTrackNumber;
            double BinEffError = sqrt(BinEffMean * (100. - BinEffMean) / simTrackNumber);
            EfficiencyEtamaxPurity->SetBinContent(i, BinEffMean);
            //EfficiencyEtamaxPurity->SetBinError(i, BinEffError);
        }
    }

    // get Eff2simPt
    for(int PtIndex = 1; PtIndex <= (int)(PtScale/2); PtIndex++) {
        double ParticleBinValue = Particle2simPtHist->GetBinContent(PtIndex);
        double STABinValue = STA2simPtHist->GetBinContent(PtIndex);
        if(ParticleBinValue == 0.)
            ParticleBinValue += 1.;
        double EfficiencyBinValue = STABinValue / ParticleBinValue * 100.;
        double EfficiencyBinError = sqrt(EfficiencyBinValue * (100. - EfficiencyBinValue) / ParticleBinValue);
        Eff2simPtHist->SetBinContent(PtIndex, EfficiencyBinValue);
        Eff2simPtHist->SetBinError(PtIndex, EfficiencyBinError);
    }
    
    double minX = 0;
    double minY = 0;
    double maxX = 110;
    double maxY = 40;

    TCanvas* OutputCanvas = new TCanvas("Canvas", "Canvas", 800, 600);
    OutputCanvas->cd();
    TPad* myPad = new TPad("Pad", "Pad", 0, 0, 1, 1);
    myPad->Draw();
    myPad->cd();

    string SaveName;

    EfficiencyforTrackingParticle->Draw("");
    SaveName=FileName+"_EfficiencyforTrackingParticle.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(1);
    maxPurityperTrackingParticle->Draw("");
    SaveName=FileName+"_maxPurityperTrackingParticle.png";OutputCanvas->SaveAs(SaveName.c_str());
    MultiplicityperTrackingPartile->Draw("");
    SaveName=FileName+"_MultiplicityperTrackingPartile.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(0);
    recTrackrefMomentummaxPurity->Draw("");
    SaveName=FileName+"_recTrackrefMomentummaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackrefPhimaxPurity->Draw("");
    SaveName=FileName+"_recTrackrefPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackrefEtamaxPurity->Draw("");
    SaveName=FileName+"_recTrackrefEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerMomentummaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerMomentummaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerPhimaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerEtamaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(1);
    recTrackinnerValidmaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerValidmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(0);
    recTrackouterMomentummaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterMomentummaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackouterPhimaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackouterEtamaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(1);
    recTrackouterValidmaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterValidmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(0);
    simTrackMomentumPtmaxPurity->Draw("");
    SaveName=FileName+"_simTrackMomentumPtmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackPhimaxPurity->Draw("");
    SaveName=FileName+"_simTrackPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackEtamaxPurity->Draw("");
    SaveName=FileName+"_simTrackEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackinnerMomentummaxPurity->Draw("");
    SaveName=FileName+"_simTrackinnerMomentummaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackinnerPhimaxPurity->Draw("");
    SaveName=FileName+"_simTrackinnerPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackinnerEtamaxPurity->Draw("");
    SaveName=FileName+"_simTrackinnerEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackinnerMatchmaxPurity->Draw("");
    SaveName=FileName+"_simTrackinnerMatchmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackouterMomentummaxPurity->Draw("");
    SaveName=FileName+"_simTrackouterMomentummaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackouterPhimaxPurity->Draw("");
    SaveName=FileName+"_simTrackouterPhimaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackouterEtamaxPurity->Draw("");
    SaveName=FileName+"_simTrackouterEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    simTrackouterMatchmaxPurity->Draw("");
    SaveName=FileName+"_simTrackouterMatchmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerMomentumofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerMomentumofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerPhiofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerPhiofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackinnerEtaofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerEtaofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(1);
    recTrackinnerValidofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackinnerValidofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(0);
    recTrackouterMomentumofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterMomentumofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackouterPhiofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterPhiofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackouterEtaofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterEtaofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(1);
    recTrackouterValidofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackouterValidofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    myPad->SetLogy(0);
    recTrackimpactMomentumofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackimpactMomentumofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackimpactPhiofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackimpactPhiofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackimpactEtaofTSOSmaxPurity->Scale(100./(double)recTrackimpactEtaofTSOSmaxPurity->GetEntries());
    recTrackimpactEtaofTSOSmaxPurity->GetXaxis()->SetTitle("Eta");
    recTrackimpactEtaofTSOSmaxPurity->GetXaxis()->CenterTitle(1);
    recTrackimpactEtaofTSOSmaxPurity->GetYaxis()->SetTitle("Fraction %");
    recTrackimpactEtaofTSOSmaxPurity->GetYaxis()->CenterTitle(1);
    recTrackimpactEtaofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackimpactEtaofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    recTrackimpactValidofTSOSmaxPurity->Draw("");
    SaveName=FileName+"_recTrackimpactValidofTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    deltaPtatimpactTSOSmaxPurity->Draw("");
    SaveName=FileName+"_deltaPtatimpactTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    Particle2simPtHist->Draw("");
    SaveName=FileName+"_Particle2simPt.png";OutputCanvas->SaveAs(SaveName.c_str());
    STA2simPtHist->Draw("");
    SaveName=FileName+"_STA2simPt.png";OutputCanvas->SaveAs(SaveName.c_str());
    Eff2simPtHist->GetXaxis()->SetTitle("simPt /Gev");
    Eff2simPtHist->GetXaxis()->CenterTitle(1);
    Eff2simPtHist->GetYaxis()->SetTitle("Efficiency %");
    Eff2simPtHist->GetYaxis()->CenterTitle(1);
    Eff2simPtHist->Draw("");
    SaveName=FileName+"_Eff2simPt.png";OutputCanvas->SaveAs(SaveName.c_str());
    EfficiencyEtamaxPurity->GetXaxis()->SetTitle("Eta");
    EfficiencyEtamaxPurity->GetXaxis()->CenterTitle(1);
    EfficiencyEtamaxPurity->GetYaxis()->SetTitle("Efficiency %");
    EfficiencyEtamaxPurity->GetYaxis()->CenterTitle(1);
    EfficiencyEtamaxPurity->Draw();
    SaveName=FileName+"_EfficiencyEtamaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
    ChargeCheckHistTSOSmaxPurityHist->GetXaxis()->SetTitle("simCharge * recCharge");
    ChargeCheckHistTSOSmaxPurityHist->GetXaxis()->CenterTitle(1);
    ChargeCheckHistTSOSmaxPurityHist->GetYaxis()->SetTitle("Fraction %");
    ChargeCheckHistTSOSmaxPurityHist->GetYaxis()->CenterTitle(1);
    ChargeCheckHistTSOSmaxPurityHist->Scale(100./double(ChargeCheckHistTSOSmaxPurityHist->GetEntries()));
    ChargeCheckHistTSOSmaxPurityHist->Draw();
    SaveName=FileName+"_ChargeCheckHistTSOSmaxPurity.png";OutputCanvas->SaveAs(SaveName.c_str());
}
