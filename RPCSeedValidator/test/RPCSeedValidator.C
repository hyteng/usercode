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

    if(debug) cout << FileName << endl;
    string theFileName = FileName;
    TFile* RootFile = TFile::Open(theFileName.c_str());
    string OutputPlotNameFix = ".png";
    string FinalOutput = "Filter_" + FileName + "_";

    unsigned int EventNumber;
    int SimTrackId;
    int SimTrackType;
    double SimTrackMomentum;
    double SimTrackDirectionPhi;
    int SimTrackCharge;
    int SimTrackValid;
    bool PassSegmentFilter;
    double SimMomentumatRef;
    double SimDirectionPhiatRef;
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
    T0->SetBranchAddress("SimTrackValid", &SimTrackValid);
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

    TH1D* SimTrackValidHist = (TH1D*) new TH1D("SimTrackValidHist", "SimTrackValidHist", 2, 0, 2);
    TH1D* SeedPtforSimTrackValidHist = (TH1D*) new TH1D("SeedPtforSimTrackValidHist", "SeedPtforSimTrackValidHist", 1000, 0, 100);
    TH1D* SeeddeltaPtforSimTrackValidHist = (TH1D*) new TH1D("SeeddeltaPtforSimTrackValidHist", "SeeddeltaPtforSimTrackValidHist", 400, -20, 20);
    TH1D* SeedPurityforSimTrackValidHist = (TH1D*) new TH1D("SeedPurityforSimTrackValidHist", "SeedPurityforSimTrackValidHist", 20, 0, 2);
    TH1D* ChargeCheckforSimTrackValidHist = (TH1D*) new TH1D("ChargeCheckforSimTrackValidHist", "ChargeCheckforSimTrackValidHist", 5, -2.5, 2.5);
    TH1D* SeedNumberforSimTrackValidHist = (TH1D*) new TH1D("SeedNumberforSimTrackValidHist", "SeedNumberforSimTrackValidHist", 30, 0, 30);
    TH1D* SeedEfficiencyforSimTrackValidHist = (TH1D*) new TH1D("SeedEfficiencyforSimTrackValidHist", "SeedEfficiencyforSimTrackValidHist", 2, 0, 2);
    TH1D* SeedNumberforSimTrackinValidHist = (TH1D*) new TH1D("SeedNumberforSimTrackinValidHist", "SeedNumberforSimTrackinValidHist", 20, 0, 20);
    TH1D* SeedEfficiencyforSimTrackinValidHist = (TH1D*) new TH1D("SeedEfficiencyforSimTrackinValidHist", "SeedEfficiencyforSimTrackinValidHist", 2, 0, 2);
    TH1D* RecBendingLastPhiHist = (TH1D*) new TH1D("RecBendingLastPhiHist", "RecBendingLastPhiHist", 628, -3.14/2, 3.14/2);

    TH2D* RecBendingPhi2PtHist = new TH2D("RecBendingPhi2PtHist", "RecBendingPhi2PtHist", 2000, -100, 100, 628, -3.14/4, 3.14/4);
    TH2D* PtRatoofRecBendingPhiHist = new TH2D("", "", 628, -3.14/4, 3.14/4, 2000, -100, 100);
    TObjArray* SimReverseBending = (TObjArray*) new TObjArray();

    unsigned int lastSeedNumber = -1;
    unsigned int lastSimTrackValid = 0;
    bool lastPassSegmentFilter = false;
    int Nentries = T0->GetEntries();
    for(int i = 0; i < Nentries; i++) {
        T0->GetEntry(i);

        if(debug) cout << "SimTrackId: " << SimTrackId << ", SimTrackType: " << SimTrackType << ", PassSegmentFilter: " << PassSegmentFilter << ", SeedNumber: " << SeedNumber << ", SeedPurity: " << SeedPurity << ", SeedCharge: " << SeedCharge << ", SimTrackCharge: " << SimTrackCharge << ", SimBendingPhi: " << SimBendingPhi << ", RecBendingPhi: " << RecBendingPhi << ", RecBendingLastPhi: " << RecBendingLastPhi << ", RecMomentumatRef: " << RecMomentumatRef << ", SimMomentumatRef: " << SimMomentumatRef << ", lastSeedNumber: " << lastSeedNumber << endl;

        if(SeedNumber != 0) {
            if(SimTrackValid == 1 && PassSegmentFilter == true) {
                //SeedPtforSimTrackValidHist->Fill(RecMomentumatRef);
                //SeeddeltaPtforSimTrackValidHist->Fill(RecMomentumatRef-SimMomentumatRef);
                SeedPurityforSimTrackValidHist->Fill(SeedPurity);
                //ChargeCheckforSimTrackValidHist->Fill(SeedCharge*SimTrackCharge);
                //RecBendingLastPhiHist->Fill(RecBendingLastPhi*(double)SimTrackCharge);
                if(SeedPurity == 1) {
                    ChargeCheckforSimTrackValidHist->Fill(SeedCharge*SimTrackCharge);
                    SeedPtforSimTrackValidHist->Fill(RecMomentumatRef);
                    SeeddeltaPtforSimTrackValidHist->Fill(RecMomentumatRef-SimMomentumatRef);
                    RecBendingPhi2PtHist->Fill(SimMomentumatRef*SimTrackCharge, RecBendingPhi);
                    double PtRato = SimMomentumatRef / RecMomentumatRef;
                    if(debug) cout << "PtRato: " << PtRato << ", at RecBendingPhi: " << RecBendingPhi << endl;
                    PtRatoofRecBendingPhiHist->Fill(RecBendingPhi, PtRato);
                    RecBendingLastPhiHist->Fill((RecBendingLastPhi == 0. ? 0 : RecBendingLastPhi/fabs(RecBendingLastPhi))*SimTrackCharge);
                }
                if(SeedCharge*SimTrackCharge == -1) {
                    if(debug) cout << "R1: " << sqrt(SimBendingEntryPositionX*SimBendingEntryPositionX+SimBendingEntryPositionY*SimBendingEntryPositionY) << ", R2: " << sqrt(SimBendingLeavePositionX*SimBendingLeavePositionX+SimBendingLeavePositionY*SimBendingLeavePositionY) << endl;
                    TLine* SimSegment = new TLine(SimBendingEntryPositionX, SimBendingEntryPositionY, SimBendingLeavePositionX, SimBendingLeavePositionY);
                    SimReverseBending->AddLast(SimSegment);
                    //RecBendingLastPhiHist->Fill(RecBendingLastPhi*SimTrackCharge);
                }
            }
        }
        else {
            if(lastSeedNumber != -1) {
                if(lastSimTrackValid == 1 && lastPassSegmentFilter == true) {
                    if(debug) cout << "Filling Valid track efficiency " << lastSeedNumber << endl;
                    SeedNumberforSimTrackValidHist->Fill(lastSeedNumber);
                    SeedEfficiencyforSimTrackValidHist->Fill(lastSeedNumber>0?1:0);
                }
                else {
                    SeedNumberforSimTrackinValidHist->Fill(lastSeedNumber);
                    SeedEfficiencyforSimTrackinValidHist->Fill(lastSeedNumber>0?1:0);
                }
            }
        }
        lastPassSegmentFilter = PassSegmentFilter;
        lastSeedNumber = SeedNumber;
        lastSimTrackValid = SimTrackValid;
    }
    if(lastSeedNumber != -1) {
        if(lastSimTrackValid == 1 && lastPassSegmentFilter == true) {
            SeedNumberforSimTrackValidHist->Fill(lastSeedNumber);
            SeedEfficiencyforSimTrackValidHist->Fill(lastSeedNumber>0?1:0);
        }
        else {
            SeedNumberforSimTrackinValidHist->Fill(lastSeedNumber);
            SeedEfficiencyforSimTrackinValidHist->Fill(lastSeedNumber>0?1:0);
        }
    }

    TCanvas* SimTrackValidCanvas = new TCanvas("SimTrackValidCanvas", "SimTrackValidCanvas", 800, 600);
    SimTrackValidCanvas->cd();
    SimTrackValidHist->Draw();
    string SimTrackValidCanvasName = FinalOutput + "SimTrackValid" + OutputPlotNameFix;
    SimTrackValidCanvas->SaveAs(SimTrackValidCanvasName.c_str());

    TCanvas* SeedPtforSimTrackValidCanvas = new TCanvas("SeedPtforSimTrackValidCanvas", "SeedPtforSimTrackValidCanvas", 800, 600);
    SeedPtforSimTrackValidCanvas->cd();
    SeedPtforSimTrackValidHist->Draw();
    string SeedPtforSimTrackValidCanvasName = FinalOutput + "SeedPtforSimTrackValid" + OutputPlotNameFix;
    SeedPtforSimTrackValidCanvas->SaveAs(SeedPtforSimTrackValidCanvasName.c_str());

    TCanvas* SeeddeltaPtforSimTrackValidCanvas = new TCanvas("SeeddeltaPtforSimTrackValidCanvas", "SeeddeltaPtforSimTrackValidCanvas", 800, 600);
    SeeddeltaPtforSimTrackValidCanvas->cd();
    SeeddeltaPtforSimTrackValidHist->Draw();
    string SeeddeltaPtforSimTrackValidCanvasName = FinalOutput + "SeeddeltaPtforSimTrackValid" + OutputPlotNameFix;
    SeeddeltaPtforSimTrackValidCanvas->SaveAs(SeeddeltaPtforSimTrackValidCanvasName.c_str());

    TCanvas* SeedPurityforSimTrackValidCanvas = new TCanvas("SeedPurityforSimTrackValidCanvas", "SeedPurityforSimTrackValidCanvas", 800, 600);
    SeedPurityforSimTrackValidCanvas->cd();
    SeedPurityforSimTrackValidHist->Draw();
    string SeedPurityforSimTrackValidCanvasName = FinalOutput + "SeedPurityforSimTrackValid" + OutputPlotNameFix;
    SeedPurityforSimTrackValidCanvas->SaveAs(SeedPurityforSimTrackValidCanvasName.c_str());

    TCanvas* ChargeCheckforSimTrackValidCanvas = new TCanvas("ChargeCheckforSimTrackValidCanvas", "ChargeCheckforSimTrackValidCanvas", 800, 600);
    ChargeCheckforSimTrackValidCanvas->cd();
    ChargeCheckforSimTrackValidHist->Draw();
    string ChargeCheckforSimTrackValidCanvasName = FinalOutput + "ChargeCheckforSimTrackValid" + OutputPlotNameFix;
    ChargeCheckforSimTrackValidCanvas->SaveAs(ChargeCheckforSimTrackValidCanvasName.c_str());

    TCanvas* SeedNumberforSimTrackValidCanvas = new TCanvas("SeedNumberforSimTrackValidCanvas", "SeedNumberforSimTrackValidCanvas", 800, 600);
    SeedNumberforSimTrackValidCanvas->cd();
    SeedNumberforSimTrackValidHist->Draw();
    string SeedNumberforSimTrackValidCanvasName = FinalOutput + "SeedNumberforSimTrackValid" + OutputPlotNameFix;
    SeedNumberforSimTrackValidCanvas->SaveAs(SeedNumberforSimTrackValidCanvasName.c_str());

    TCanvas* SeedEfficiencyforSimTrackValidCanvas = new TCanvas("SeedEfficiencyforSimTrackValidCanvas", "SeedEfficiencyforSimTrackValidCanvas", 800, 600);
    SeedEfficiencyforSimTrackValidCanvas->cd();
    SeedEfficiencyforSimTrackValidHist->Draw();
    string SeedEfficiencyforSimTrackValidCanvasName = FinalOutput + "SeedEfficiencyforSimTrackValid" + OutputPlotNameFix;
    SeedEfficiencyforSimTrackValidCanvas->SaveAs(SeedEfficiencyforSimTrackValidCanvasName.c_str());

    TCanvas* SeedNumberforSimTrackinValidCanvas = new TCanvas("SeedNumberforSimTrackinValidCanvas", "SeedNumberforSimTrackinValidCanvas", 800, 600);
    SeedNumberforSimTrackinValidCanvas->cd();
    SeedNumberforSimTrackinValidHist->Draw();
    string SeedNumberforSimTrackinValidCanvasName = FinalOutput + "SeedNumberforSimTrackinValid" + OutputPlotNameFix;
    SeedNumberforSimTrackinValidCanvas->SaveAs(SeedNumberforSimTrackinValidCanvasName.c_str());

    TCanvas* SeedEfficiencyforSimTrackinValidCanvas = new TCanvas("SeedEfficiencyforSimTrackinValidCanvas", "SeedEfficiencyforSimTrackinValidCanvas", 800, 600);
    SeedEfficiencyforSimTrackinValidCanvas->cd();
    SeedEfficiencyforSimTrackinValidHist->Draw();
    string SeedEfficiencyforSimTrackinValidCanvasName = FinalOutput + "SeedEfficiencyforSimTrackinValid" + OutputPlotNameFix;
    SeedEfficiencyforSimTrackinValidCanvas->SaveAs(SeedEfficiencyforSimTrackinValidCanvasName.c_str());

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
