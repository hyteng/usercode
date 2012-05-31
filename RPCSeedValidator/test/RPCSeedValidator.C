#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TMath.h"

#define debug 1

void RPCSeedValidator(string FileName) {

    gStyle->SetOptStat("");
    gStyle->SetOptTitle(0);
    if(debug) cout << FileName << endl;
    string theFileName = "../" + FileName;
    TFile* RootFile = TFile::Open(theFileName.c_str());
    string OutputPlotNameFix = ".eps";
    string FinalOutput = FileName + "_";

    unsigned int EventNumber;
    int SimTrackId;
    int SimTrackType;
    double SimTrackMomentum;
    double SimTrackDirectionPhi;
    int SimTrackCharge;
    int SimTrackvalid;
    bool PassSegmentFilter;
    double SimMomentumatRef;
    double SimDirectionPhiatRef;
    double SimDirectionEtaatRef;
    double SimBendingPhi;
    double SimBendingEntryPositionX;
    double SimBendingEntryPositionY;
    double SimBendingEntryPositionZ;
    double SimBendingLeavePositionX;
    double SimBendingLeavePositionY;
    double SimBendingLeavePositionZ;
    unsigned int SeedNumber;
    int SeedCharge;
    double SeedPurity;
    double SeedQuality;
    double RecMomentumatRef;
    double RecDirectionPhiatRef;
    double RecDirectionEtaatRef;
    double RecBendingPhi;
    double RecBendingEntryPositionX;
    double RecBendingEntryPositionY;
    double RecBendingEntryPositionZ;
    double RecBendingLeavePositionX;
    double RecBendingLeavePositionY;
    double RecBendingLeavePositionZ;
    double RecBendingLastPhi;

    TTree* T0 = (TTree*)RootFile->Get("ExTree");
    T0->SetBranchAddress("EventNumber", &EventNumber);
    T0->SetBranchAddress("SimTrackId", &SimTrackId);
    T0->SetBranchAddress("SimTrackType", &SimTrackType);
    T0->SetBranchAddress("SimTrackMomentum", &SimTrackMomentum);
    T0->SetBranchAddress("SimTrackDirectionPhi", &SimTrackDirectionPhi);
    T0->SetBranchAddress("SimTrackCharge", &SimTrackCharge);
    T0->SetBranchAddress("SimTrackValid", &SimTrackvalid);
    T0->SetBranchAddress("PassSegmentFilter", &PassSegmentFilter);
    T0->SetBranchAddress("SimMomentumatRef", &SimMomentumatRef);
    T0->SetBranchAddress("SimDirectionPhiatRef", &SimDirectionPhiatRef);
    T0->SetBranchAddress("SimBendingPhi", &SimBendingPhi);
    T0->SetBranchAddress("SimBendingEntryPositionX", &SimBendingEntryPositionX);
    T0->SetBranchAddress("SimBendingEntryPositionY", &SimBendingEntryPositionY);
    T0->SetBranchAddress("SimBendingEntryPositionZ", &SimBendingEntryPositionZ);
    T0->SetBranchAddress("SimBendingLeavePositionX", &SimBendingLeavePositionX);
    T0->SetBranchAddress("SimBendingLeavePositionY", &SimBendingLeavePositionY);
    T0->SetBranchAddress("SimBendingLeavePositionZ", &SimBendingLeavePositionZ);
    T0->SetBranchAddress("SeedNumber", &SeedNumber);
    T0->SetBranchAddress("SeedCharge", &SeedCharge);
    T0->SetBranchAddress("SeedPurity", &SeedPurity);
    T0->SetBranchAddress("SeedQuality", &SeedQuality);
    T0->SetBranchAddress("RecMomentumatRef", &RecMomentumatRef);
    T0->SetBranchAddress("RecDirectionPhiatRef", &RecDirectionPhiatRef);
    T0->SetBranchAddress("RecBendingPhi", &RecBendingPhi);
    T0->SetBranchAddress("RecBendingEntryPositionX", &RecBendingEntryPositionX);
    T0->SetBranchAddress("RecBendingEntryPositionY", &RecBendingEntryPositionY);
    T0->SetBranchAddress("RecBendingEntryPositionZ", &RecBendingEntryPositionZ);
    T0->SetBranchAddress("RecBendingLeavePositionX", &RecBendingLeavePositionX);
    T0->SetBranchAddress("RecBendingLeavePositionY", &RecBendingLeavePositionY);
    T0->SetBranchAddress("RecBendingLeavePositionZ", &RecBendingLeavePositionZ);
    T0->SetBranchAddress("RecBendingLastPhi", &RecBendingLastPhi);

    TH1D* SimTrackvalidHist = (TH1D*) new TH1D("SimTrackvalidHist", "SimTrackvalidHist", 2, 0, 2);
    TH1D* SeedPtforSimTrackvalidHist = (TH1D*) new TH1D("SeedPtforSimTrackvalidHist", "SeedPtforSimTrackvalidHist", 100, 0, 100);
    TH1D* SeeddeltaPtforSimTrackvalidHist = (TH1D*) new TH1D("SeeddeltaPtforSimTrackvalidHist", "SeeddeltaPtforSimTrackvalidHist", 150, -3., 3.);
    TH1D* SeeddeltaPhiforSimTrackvalidHist = (TH1D*) new TH1D("SeeddeltaPhiforSimTrackvalidHist", "SeeddeltaPhiforSimTrackvalidHist", 628, -3.14/6, 3.14/6);
    TH1D* SeeddeltaEtaforSimTrackvalidHist = (TH1D*) new TH1D("SeeddeltaEtaforSimTrackvalidHist", "SeeddeltaEtaforSimTrackvalidHist", 200, -1., 1.);
    TH1D* SeedPurityforSimTrackvalidHist = (TH1D*) new TH1D("SeedPurityforSimTrackvalidHist", "SeedPurityforSimTrackvalidHist", 20, 0, 2);
    TH1D* ChargeCheckforSimTrackvalidHist = (TH1D*) new TH1D("ChargeCheckforSimTrackvalidHist", "ChargeCheckforSimTrackvalidHist", 5, -2.5, 2.5);
    TH1D* SeedNumberforSimTrackvalidHist = (TH1D*) new TH1D("SeedNumberforSimTrackvalidHist", "SeedNumberforSimTrackvalidHist", 30, 0, 30);
    TH1D* SeedEfficiencyforSimTrackvalidHist = (TH1D*) new TH1D("SeedEfficiencyforSimTrackvalidHist", "SeedEfficiencyforSimTrackvalidHist", 2, 0, 2);
    TH1D* SeedNumberforSimTrackinvalidHist = (TH1D*) new TH1D("SeedNumberforSimTrackinvalidHist", "SeedNumberforSimTrackinvalidHist", 20, 0, 20);
    TH1D* SeedEfficiencyforSimTrackinvalidHist = (TH1D*) new TH1D("SeedEfficiencyforSimTrackinvalidHist", "SeedEfficiencyforSimTrackinvalidHist", 2, 0, 2);
    TH1D* RecBendingLastPhiHist = (TH1D*) new TH1D("RecBendingLastPhiHist", "RecBendingLastPhiHist", 628, -3.14/2, 3.14/2);
    TH1D* SeedEfficiencyHist = (TH1D*) new TH1D("SeedEfficiency", "SeedEfficiency", 2, 0, 2);
    TH2D* RecBendingPhi2PtHist = new TH2D("RecBendingPhi2PtHist", "RecBendingPhi2PtHist", 2000, -100, 100, 628, -3.14/4, 3.14/4);
    TH2D* PtRatoofRecBendingPhiHist = new TH2D("", "", 628, -3.14/4, 3.14/4, 2000, -100, 100);
    TObjArray* SimReverseBending = (TObjArray*) new TObjArray();

    unsigned int LastSeedNumber = -1;
    unsigned int LastSimTrackvalid = 0;
    bool LastPassSegmentFilter = false;
    bool LastPurityFull = false;
    int Nentries = T0->GetEntries();
    for(int i = 0; i < Nentries; i++) {
        T0->GetEntry(i);

        if(debug) cout << "SimTrackId: " << SimTrackId << ", SimTrackType: " << SimTrackType << ", PassSegmentFilter: " << PassSegmentFilter << ", SeedNumber: " << SeedNumber << ", SeedPurity: " << SeedPurity << ", SeedCharge: " << SeedCharge << ", SimTrackCharge: " << SimTrackCharge << ", SimBendingPhi: " << SimBendingPhi << ", RecBendingPhi: " << RecBendingPhi << ", RecBendingLastPhi: " << RecBendingLastPhi << ", RecMomentumatRef: " << RecMomentumatRef << ", SimMomentumatRef: " << SimMomentumatRef << ", LastSeedNumber: " << LastSeedNumber << endl;

        //if(SimTrackMomentum > 20.0)
            //continue;

        if(SeedNumber != 0) {
            if(SimTrackvalid == 1 && SeedPurity == 1. && PassSegmentFilter == true) {
                SeedPtforSimTrackvalidHist->Fill(RecMomentumatRef);
                SeeddeltaPtforSimTrackvalidHist->Fill((RecMomentumatRef-SimMomentumatRef)/SimMomentumatRef);
                SeeddeltaPhiforSimTrackvalidHist->Fill(RecDirectionPhiatRef-SimDirectionPhiatRef);
                SeeddeltaEtaforSimTrackvalidHist->Fill(RecDirectionEtaatRef-SimDirectionEtaatRef);
                RecBendingPhi2PtHist->Fill(SimMomentumatRef*SimTrackCharge, RecBendingPhi);
                double PtRato = SimMomentumatRef / RecMomentumatRef;
                if(debug) cout << "PtRato: " << PtRato << ", at RecBendingPhi: " << RecBendingPhi << endl;
                PtRatoofRecBendingPhiHist->Fill(RecBendingPhi, PtRato);
                RecBendingLastPhiHist->Fill((RecBendingLastPhi == 0. ? 0 : RecBendingLastPhi/fabs(RecBendingLastPhi))*SimTrackCharge);
                SeedPurityforSimTrackvalidHist->Fill(SeedPurity);
                ChargeCheckforSimTrackvalidHist->Fill(SeedCharge*SimTrackCharge);
                /*
                if(SeedPurity == 1) {
                    //SeedPtforSimTrackvalidHist->Fill(RecMomentumatRef);
                    //SeeddeltaPtforSimTrackvalidHist->Fill(RecMomentumatRef-SimMomentumatRef);
                    RecBendingPhi2PtHist->Fill(SimMomentumatRef*SimTrackCharge, RecBendingPhi);
                    double PtRato = SimMomentumatRef / RecMomentumatRef;
                    if(debug) cout << "PtRato: " << PtRato << ", at RecBendingPhi: " << RecBendingPhi << endl;
                    PtRatoofRecBendingPhiHist->Fill(RecBendingPhi, PtRato);
                    RecBendingLastPhiHist->Fill((RecBendingLastPhi == 0. ? 0 : RecBendingLastPhi/fabs(RecBendingLastPhi))*SimTrackCharge);
                }
                */
                if(SeedCharge*SimTrackCharge == -1) {
                    if(debug) cout << "R1: " << sqrt(SimBendingEntryPositionX*SimBendingEntryPositionX+SimBendingEntryPositionY*SimBendingEntryPositionY) << ", R2: " << sqrt(SimBendingLeavePositionX*SimBendingLeavePositionX+SimBendingLeavePositionY*SimBendingLeavePositionY) << endl;
                    TLine* SimSegment = new TLine(SimBendingEntryPositionX, SimBendingEntryPositionY, SimBendingLeavePositionX, SimBendingLeavePositionY);
                    SimReverseBending->AddLast(SimSegment);
                }
                LastPurityFull = true;
            }
        }
        else {
            if(LastSeedNumber != -1) {
                if(LastSimTrackvalid == 1 && LastPassSegmentFilter == true) {
                    if(debug) cout << "Filling valid track efficiency " << LastSeedNumber << endl;
                    SeedNumberforSimTrackvalidHist->Fill(LastSeedNumber);
                    SeedEfficiencyforSimTrackvalidHist->Fill(LastPurityFull==true?1:0);
                }
                else {
                    SeedNumberforSimTrackinvalidHist->Fill(LastSeedNumber);
                    SeedEfficiencyforSimTrackinvalidHist->Fill(LastPurityFull==true?1:0);
                }
                LastPurityFull = false;
            }
        }
        LastPassSegmentFilter = PassSegmentFilter;
        LastSeedNumber = SeedNumber;
        LastSimTrackvalid = SimTrackvalid;
    }
    if(LastSeedNumber != -1) {
        if(LastSimTrackvalid == 1 && LastPassSegmentFilter == true) {
            SeedNumberforSimTrackvalidHist->Fill(LastSeedNumber);
            SeedEfficiencyforSimTrackvalidHist->Fill(LastPurityFull==true?1:0);
        }
        else {
            SeedNumberforSimTrackinvalidHist->Fill(LastSeedNumber);
            SeedEfficiencyforSimTrackinvalidHist->Fill(LastPurityFull == true?1:0);
        }
    }

    TCanvas* SimTrackvalidCanvas = new TCanvas("SimTrackvalidCanvas", "SimTrackvalidCanvas", 800, 600);
    SimTrackvalidCanvas->cd();
    SimTrackvalidHist->Draw();
    string SimTrackvalidCanvasName = FinalOutput + "SimTrackvalid" + OutputPlotNameFix;
    SimTrackvalidCanvas->SaveAs(SimTrackvalidCanvasName.c_str());

    TCanvas* SeedPtforSimTrackvalidCanvas = new TCanvas("SeedPtforSimTrackvalidCanvas", "SeedPtforSimTrackvalidCanvas", 800, 600);
    SeedPtforSimTrackvalidCanvas->cd();
    SeedPtforSimTrackvalidHist->Draw();
    string SeedPtforSimTrackvalidCanvasName = FinalOutput + "SeedPtforSimTrackvalid" + OutputPlotNameFix;
    SeedPtforSimTrackvalidCanvas->SaveAs(SeedPtforSimTrackvalidCanvasName.c_str());

    TCanvas* SeeddeltaPtforSimTrackvalidCanvas = new TCanvas("SeeddeltaPtforSimTrackvalidCanvas", "SeeddeltaPtforSimTrackvalidCanvas", 800, 600);
    SeeddeltaPtforSimTrackvalidCanvas->cd();
    SeeddeltaPtforSimTrackvalidHist->SetStats(1);
    gStyle->SetOptFit(0111);
    //SeeddeltaPtforSimTrackvalidHist->Fit("gaus", "", "", -1., 1.);
    SeeddeltaPtforSimTrackvalidHist->GetXaxis()->SetTitle("(recPt-simPt)/simPt");
    SeeddeltaPtforSimTrackvalidHist->GetXaxis()->CenterTitle();
    SeeddeltaPtforSimTrackvalidHist->Draw();
    string SeeddeltaPtforSimTrackvalidCanvasName = FinalOutput + "SeeddeltaPtforSimTrackvalid" + OutputPlotNameFix;
    SeeddeltaPtforSimTrackvalidCanvas->SaveAs(SeeddeltaPtforSimTrackvalidCanvasName.c_str());

    TCanvas* SeeddeltaPhiforSimTrackvalidCanvas = new TCanvas("SeeddeltaPhiforSimTrackvalidCanvas", "SeeddeltaPhiforSimTrackvalidCanvas", 800, 600);
    SeeddeltaPhiforSimTrackvalidCanvas->cd();
    SeeddeltaPhiforSimTrackvalidHist->SetStats(1);
    SeeddeltaPhiforSimTrackvalidHist->GetXaxis()->SetTitle("(recPhi-simPhi)/simPhi");
    SeeddeltaPhiforSimTrackvalidHist->GetXaxis()->CenterTitle();
    SeeddeltaPhiforSimTrackvalidHist->Draw();
    string SeeddeltaPhiforSimTrackvalidCanvasName = FinalOutput + "SeeddeltaPhiforSimTrackvalid" + OutputPlotNameFix;
    SeeddeltaPhiforSimTrackvalidCanvas->SaveAs(SeeddeltaPhiforSimTrackvalidCanvasName.c_str());

    TCanvas* SeeddeltaEtaforSimTrackvalidCanvas = new TCanvas("SeeddeltaEtaforSimTrackvalidCanvas", "SeeddeltaEtaforSimTrackvalidCanvas", 800, 600);
    SeeddeltaEtaforSimTrackvalidCanvas->cd();
    SeeddeltaEtaforSimTrackvalidHist->SetStats(1);
    SeeddeltaEtaforSimTrackvalidHist->GetXaxis()->SetTitle("(recEta-simEta)/simEta");
    SeeddeltaEtaforSimTrackvalidHist->GetXaxis()->CenterTitle();
    SeeddeltaEtaforSimTrackvalidHist->Draw();
    string SeeddeltaEtaforSimTrackvalidCanvasName = FinalOutput + "SeeddeltaEtaforSimTrackvalid" + OutputPlotNameFix;
    SeeddeltaEtaforSimTrackvalidCanvas->SaveAs(SeeddeltaEtaforSimTrackvalidCanvasName.c_str());


    TCanvas* SeedPurityforSimTrackvalidCanvas = new TCanvas("SeedPurityforSimTrackvalidCanvas", "SeedPurityforSimTrackvalidCanvas", 800, 600);
    SeedPurityforSimTrackvalidCanvas->cd();
    SeedPurityforSimTrackvalidHist->Draw();
    string SeedPurityforSimTrackvalidCanvasName = FinalOutput + "SeedPurityforSimTrackvalid" + OutputPlotNameFix;
    SeedPurityforSimTrackvalidCanvas->SaveAs(SeedPurityforSimTrackvalidCanvasName.c_str());

    TCanvas* ChargeCheckforSimTrackvalidCanvas = new TCanvas("ChargeCheckforSimTrackvalidCanvas", "ChargeCheckforSimTrackvalidCanvas", 800, 600);
    ChargeCheckforSimTrackvalidCanvas->cd();
    double HistEntries = ChargeCheckforSimTrackvalidHist->GetEntries() / 100.;
    ChargeCheckforSimTrackvalidHist->Scale(1./HistEntries);
    ChargeCheckforSimTrackvalidHist->GetXaxis()->SetTitle("simCharge*recCharge");
    ChargeCheckforSimTrackvalidHist->GetXaxis()->CenterTitle(1);
    ChargeCheckforSimTrackvalidHist->GetYaxis()->SetTitle("fraction %");
    ChargeCheckforSimTrackvalidHist->GetYaxis()->CenterTitle(1);
    ChargeCheckforSimTrackvalidHist->Draw();
    string ChargeCheckforSimTrackvalidCanvasName = FinalOutput + "ChargeCheckforSimTrackvalid" + OutputPlotNameFix;
    ChargeCheckforSimTrackvalidCanvas->SaveAs(ChargeCheckforSimTrackvalidCanvasName.c_str());

    TCanvas* SeedNumberforSimTrackvalidCanvas = new TCanvas("SeedNumberforSimTrackvalidCanvas", "SeedNumberforSimTrackvalidCanvas", 800, 600);
    SeedNumberforSimTrackvalidCanvas->cd();
    SeedNumberforSimTrackvalidHist->Draw();
    string SeedNumberforSimTrackvalidCanvasName = FinalOutput + "SeedNumberforSimTrackvalid" + OutputPlotNameFix;
    SeedNumberforSimTrackvalidCanvas->SaveAs(SeedNumberforSimTrackvalidCanvasName.c_str());

    TCanvas* SeedEfficiencyforSimTrackvalidCanvas = new TCanvas("SeedEfficiencyforSimTrackvalidCanvas", "SeedEfficiencyforSimTrackvalidCanvas", 800, 600);
    SeedEfficiencyforSimTrackvalidCanvas->cd();
    SeedEfficiencyforSimTrackvalidHist->Draw();
    string SeedEfficiencyforSimTrackvalidCanvasName = FinalOutput + "SeedEfficiencyforSimTrackvalid" + OutputPlotNameFix;
    SeedEfficiencyforSimTrackvalidCanvas->SaveAs(SeedEfficiencyforSimTrackvalidCanvasName.c_str());

    TCanvas* SeedNumberforSimTrackinvalidCanvas = new TCanvas("SeedNumberforSimTrackinvalidCanvas", "SeedNumberforSimTrackinvalidCanvas", 800, 600);
    SeedNumberforSimTrackinvalidCanvas->cd();
    SeedNumberforSimTrackinvalidHist->Draw();
    string SeedNumberforSimTrackinvalidCanvasName = FinalOutput + "SeedNumberforSimTrackinvalid" + OutputPlotNameFix;
    SeedNumberforSimTrackinvalidCanvas->SaveAs(SeedNumberforSimTrackinvalidCanvasName.c_str());

    TCanvas* SeedEfficiencyforSimTrackinvalidCanvas = new TCanvas("SeedEfficiencyforSimTrackinvalidCanvas", "SeedEfficiencyforSimTrackinvalidCanvas", 800, 600);
    SeedEfficiencyforSimTrackinvalidCanvas->cd();
    SeedEfficiencyforSimTrackinvalidHist->Draw();
    string SeedEfficiencyforSimTrackinvalidCanvasName = FinalOutput + "SeedEfficiencyforSimTrackinvalid" + OutputPlotNameFix;
    SeedEfficiencyforSimTrackinvalidCanvas->SaveAs(SeedEfficiencyforSimTrackinvalidCanvasName.c_str());

    double SeedEfficiencyforSimTrackinvalid = 100. * SeedEfficiencyforSimTrackinvalidHist->GetMean();
    double SeedEfficiencyforSimTrackvalid = 100. * SeedEfficiencyforSimTrackvalidHist->GetMean();
    SeedEfficiencyHist->SetBinContent(1, SeedEfficiencyforSimTrackinvalid);
    SeedEfficiencyHist->SetBinContent(2, SeedEfficiencyforSimTrackvalid);
    SeedEfficiencyHist->GetXaxis()->SetBinLabel(1, "for invalid simTrack");
    SeedEfficiencyHist->GetXaxis()->SetBinLabel(2, "for valid simTrack");
    TCanvas* SeedEfficiencyCanvas = new TCanvas("SeedEfficiencyCanvas", "SeedEfficiencyCanvas", 800, 600);
    SeedEfficiencyCanvas->cd();
    SeedEfficiencyHist->GetYaxis()->SetTitle("Efficiency %");
    SeedEfficiencyHist->GetYaxis()->CenterTitle(1);
    SeedEfficiencyHist->SetMarkerStyle(3);
    SeedEfficiencyHist->SetMarkerSize(3);
    SeedEfficiencyHist->Draw("P");
    string SeedEfficiencyCanvasName = FinalOutput + "SeedEfficiency" + OutputPlotNameFix;
    SeedEfficiencyCanvas->SaveAs(SeedEfficiencyCanvasName.c_str());

    TCanvas* RecBendingPhi2PtCanvas = new TCanvas("RecBendingPhi2PtCanvas", "RecBendingPhi2PtCanvas", 800, 600);
    RecBendingPhi2PtCanvas->cd();
    RecBendingPhi2PtHist->Draw();
    string RecBendingPhi2PtCanvasName = FinalOutput + "RecBendingPhi2Pt" + OutputPlotNameFix;
    RecBendingPhi2PtCanvas->SaveAs(RecBendingPhi2PtCanvasName.c_str());

    TCanvas* PtRatoofRecBendingPhiCanvas = new TCanvas("PtRatoofRecBendingPhiCanvas", "PtRatoofRecBendingPhiCanvas", 800, 600);
    PtRatoofRecBendingPhiCanvas->cd();
    PtRatoofRecBendingPhiHist->Draw();
    string PtRatoofRecBendingPhiCanvasName = FinalOutput + "PtRatoofRecBendingPhi" + OutputPlotNameFix;
    PtRatoofRecBendingPhiCanvas->SaveAs(PtRatoofRecBendingPhiCanvasName.c_str());

    TCanvas* RecBendingLastPhiCanvas = new TCanvas("RecBendingLastPhiCanvas", "RecBendingLastPhiCanvas", 800, 600);
    RecBendingLastPhiCanvas->cd();
    RecBendingLastPhiHist->Draw();
    string RecBendingLastPhiCanvasName = FinalOutput + "RecBendingLastPhi" + OutputPlotNameFix;
    RecBendingLastPhiCanvas->SaveAs(RecBendingLastPhiCanvasName.c_str());

    Int_t linsav = gStyle->GetLineWidth();
    gStyle->SetLineWidth(2);
    TCanvas* SimReverseBendingCanvas = new TCanvas("SimReverseBendingCanvas", "SimReverseBendingCanvas", 800, 800);
    SimReverseBendingCanvas->cd();
    TPad* SimReverseBendingPad = new TPad("SimReverseBendingPad", "SimReverseBendingPad", 0, 0, 1, 1);
    SimReverseBendingPad->Draw();
    SimReverseBendingPad->cd();
    SimReverseBendingPad->Range(-800, -800, 800, 800);
    unsigned int segmentNumber = SimReverseBending->GetEntries();
    cout << "Number of segments: " << segmentNumber << endl;
    for(unsigned int j = 0; j < segmentNumber; j++) {
        ((TLine*)(SimReverseBending->At(j)))->Print();
        ((TLine*)(SimReverseBending->At(j)))->Draw("SAME");
    }
    string SimReverseBendingCanvasName = FinalOutput + "SimReverseBending" + OutputPlotNameFix;
    SimReverseBendingCanvas->SaveAs(SimReverseBendingCanvasName.c_str());
}
