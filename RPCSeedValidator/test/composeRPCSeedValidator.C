#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TMath.h"

#define debug 0

void composeRPCSeedValidator(string FileListName, int FileNumber) {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat("");


    if(debug) cout << FileListName << endl;
    string theFileName;
    ifstream composeFileList;
    composeFileList.open(FileListName.c_str());

    string OutputPlotNameFix = ".eps";

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

    vector<double> simPtArray(0);
    vector<double> EfficiencyArray(0);
    vector<double> InverseChargeRatoArray(0);
    TObjArray* myChargeCheckHist = new TObjArray();
    TObjArray* myDeltaPtHist = new TObjArray();
    TObjArray* myEffonValidHist = new TObjArray();

    for(int Index = 0; Index < FileNumber; Index++) {
        getline(composeFileList, theFileName);
        string Output = theFileName + "_";
        string FullFilePath = "hist/" + theFileName;
        TFile* RootFile = TFile::Open(FullFilePath.c_str());
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

        string SeeddeltaPtforSimTrackValidHistName = Output + "SeeddeltaPtforSimTrackValidHist";
        TH1D* SeeddeltaPtforSimTrackValidHist = (TH1D*) new TH1D(SeeddeltaPtforSimTrackValidHistName.c_str(), SeeddeltaPtforSimTrackValidHistName.c_str(), 200, -10, 10);
        string ChargeCheckforSimTrackValidHistName = Output + "ChargeCheckforSimTrackValidHist";
        TH1D* ChargeCheckforSimTrackValidHist = (TH1D*) new TH1D(ChargeCheckforSimTrackValidHistName.c_str(), ChargeCheckforSimTrackValidHistName.c_str(), 5, -2.5, 2.5);
        string SeedEfficiencyforSimTrackValidHistName = Output + "SeedEfficiencyforSimTrackValidHist";
        TH1D* SeedEfficiencyforSimTrackValidHist = (TH1D*) new TH1D(SeedEfficiencyforSimTrackValidHistName.c_str(), SeedEfficiencyforSimTrackValidHistName.c_str(), 2, 0, 2);

        unsigned int lastSeedNumber = -1;
        unsigned int lastSimTrackValid = 0;
        bool lastPassSegmentFilter = false;
        int Nentries = T0->GetEntries();
        for(int i = 0; i < Nentries; i++) {
            T0->GetEntry(i);

            if(debug) cout << "SimTrackId: " << SimTrackId << ", SimTrackType: " << SimTrackType << ", PassSegmentFilter: " << PassSegmentFilter << ", SeedNumber: " << SeedNumber << ", SeedPurity: " << SeedPurity << ", SeedCharge: " << SeedCharge << ", SimTrackCharge: " << SimTrackCharge << ", SimBendingPhi: " << SimBendingPhi << ", RecBendingPhi: " << RecBendingPhi << ", RecBendingLastPhi: " << RecBendingLastPhi << ", RecMomentumatRef: " << RecMomentumatRef << ", SimMomentumatRef: " << SimMomentumatRef << ", lastSeedNumber: " << lastSeedNumber << endl;

            if(SeedNumber != 0) {
                if(SimTrackValid == 1 && PassSegmentFilter == true) {
                    ChargeCheckforSimTrackValidHist->Fill(SeedCharge*SimTrackCharge);
                    SeeddeltaPtforSimTrackValidHist->Fill((RecMomentumatRef-SimMomentumatRef)/SimMomentumatRef);
                }
            }
            else {
                if(lastSeedNumber != -1) {
                    if(lastSimTrackValid == 1 && lastPassSegmentFilter == true) {
                        if(debug) cout << "Filling Valid track efficiency " << lastSeedNumber << endl;
                        SeedEfficiencyforSimTrackValidHist->Fill(lastSeedNumber>0?1:0);
                    }
                }
            }
            lastPassSegmentFilter = PassSegmentFilter;
            lastSeedNumber = SeedNumber;
            lastSimTrackValid = SimTrackValid;
        }
        if(lastSeedNumber != -1) {
            if(lastSimTrackValid == 1 && lastPassSegmentFilter == true) {
                SeedEfficiencyforSimTrackValidHist->Fill(lastSeedNumber>0?1:0);
            }
        }

        double EfficiencyforPt = SeedEfficiencyforSimTrackValidHist->GetMean();
        EfficiencyforPt *= 100.;
        double InverseChargeEvent = ChargeCheckforSimTrackValidHist->GetBinContent(2);
        double CoverseChargeEvent = ChargeCheckforSimTrackValidHist->GetBinContent(4);
        double InverseChargeRato = 100. * InverseChargeEvent / (CoverseChargeEvent + InverseChargeEvent);
        simPtArray.push_back(SimTrackMomentum);
        EfficiencyArray.push_back(EfficiencyforPt);
        InverseChargeRatoArray.push_back(InverseChargeRato);
        if(debug) cout << SimTrackMomentum << ", " << EfficiencyforPt << ", " << InverseChargeRato << endl;

        myChargeCheckHist->AddLast(ChargeCheckforSimTrackValidHist);
        myDeltaPtHist->AddLast(SeeddeltaPtforSimTrackValidHist);
        myEffonValidHist->AddLast(SeedEfficiencyforSimTrackValidHist);
    }

    TH1D* ScanDeltaPtforSimTrackValidHist = (TH1D*) new TH1D("ScanDeltaPtforSimTrackValid", "ScanDeltaPtforSimTrackValid", 505, 0, 101);
    TH1D* ScanChargeCheckforSimTrackValidHist = (TH1D*) new TH1D("ScanInverseChargeCheckforSimTrackValid", "ScanInverseChargeCheckforSimTrackValid", 505, 0, 101);
    TH1D* ScanEfficienyforSimTrackValidHist = (TH1D*) new TH1D("ScanEfficienyforSimTrackValid", "ScanEfficienyforSimTrackValid", 505, 0, 101);


    for(int Index = 0; Index < FileNumber; Index++) {
        double EfficiencyforPt = ((TH1D*)(myEffonValidHist->At(Index)))->GetMean();
        double InverseChargeEvent = ((TH1D*)(myChargeCheckHist->At(Index)))->GetBinContent(2);
        double CoverseChargeEvent = ((TH1D*)(myChargeCheckHist->At(Index)))->GetBinContent(4);
        double InverseChargeRato = 100. * InverseChargeEvent / (CoverseChargeEvent + InverseChargeEvent);
        if(((TH1D*)(myDeltaPtHist->At(Index)))->GetEntries() == 0)
            continue;

        //TCanvas* FitCanvas = new TCanvas("FitCanvas", "FitCanvas", 800, 600);
        //FitCanvas->cd();
        int MaxBin = ((TH1D*)(myDeltaPtHist->At(Index)))->GetMaximumBin();
        double MaxBinCenter = ((TH1D*)(myDeltaPtHist->At(Index)))->GetBinCenter(MaxBin);
        double RMS = ((TH1D*)(myDeltaPtHist->At(Index)))->GetRMS();
        ((TH1D*)(myDeltaPtHist->At(Index)))->Fit("gaus", "", "", MaxBinCenter-0.5*RMS, MaxBinCenter+0.5*RMS);
        double DeltaPtValue = ((TH1D*)(myDeltaPtHist->At(Index)))->GetFunction("gaus")->GetParameter(1);
        double DeltaPtError = ((TH1D*)(myDeltaPtHist->At(Index)))->GetFunction("gaus")->GetParameter(2); 
        if(DeltaPtValue < (MaxBinCenter-RMS)) {
            DeltaPtValue = MaxBinCenter;
            DeltaPtError = RMS;
        }

        std::stringstream TempIndex;
        TempIndex << Index;
        string FitIndex = TempIndex.str();
        string FitSaveName = "Fit" + FitIndex + OutputPlotNameFix;
        //FitCanvas->SaveAs(FitSaveName.c_str());

        int PtBinNumber = ScanDeltaPtforSimTrackValidHist->FindBin(simPtArray[Index]);
        ScanDeltaPtforSimTrackValidHist->SetBinContent(PtBinNumber, DeltaPtValue);
        ScanDeltaPtforSimTrackValidHist->SetBinError(PtBinNumber, DeltaPtError);
        ScanChargeCheckforSimTrackValidHist->SetBinContent(PtBinNumber, InverseChargeRato);
        ScanEfficienyforSimTrackValidHist->SetBinContent(PtBinNumber, EfficiencyforPt*100.);
    }

    TGraph* ScanEfficiencyforSimTrackValidGraph = (TGraph*)new TGraph(simPtArray.size(), &simPtArray[0], &EfficiencyArray[0]);
    TGraph* ScanChargeCheckforSimTrackValidGraph = (TGraph*)new TGraph(simPtArray.size(), &simPtArray[0], &InverseChargeRatoArray[0]);
    ScanEfficiencyforSimTrackValidGraph->SetMarkerColor(2);ScanEfficiencyforSimTrackValidGraph->SetLineColor(2);ScanEfficiencyforSimTrackValidGraph->SetFillColor(0);ScanEfficiencyforSimTrackValidGraph->SetMarkerStyle(3);ScanEfficiencyforSimTrackValidGraph->SetLineStyle(2);
    ScanChargeCheckforSimTrackValidGraph->SetMarkerColor(2);ScanChargeCheckforSimTrackValidGraph->SetLineColor(2);ScanChargeCheckforSimTrackValidGraph->SetFillColor(0);ScanChargeCheckforSimTrackValidGraph->SetMarkerStyle(3);ScanChargeCheckforSimTrackValidGraph->SetLineStyle(2);

    TCanvas* OutputCanvas0 = new TCanvas("Canvas0", "Canvas0", 800, 600);
    OutputCanvas0->cd();

    ScanDeltaPtforSimTrackValidHist->SetStats(0);
    //ScanDeltaPtforSimTrackValidHist->SetTitle("ScanDeltaPtforSimTrackValid");
    ScanDeltaPtforSimTrackValidHist->SetLineColor(kRed);
    ScanDeltaPtforSimTrackValidHist->GetXaxis()->SetTitle("simPt/Gev");
    ScanDeltaPtforSimTrackValidHist->GetXaxis()->CenterTitle(1);
    ScanDeltaPtforSimTrackValidHist->GetYaxis()->SetTitle("(seedPt-simPt)/simPt");
    ScanDeltaPtforSimTrackValidHist->GetYaxis()->CenterTitle(1);
    ScanDeltaPtforSimTrackValidHist->Draw("");
    string SaveName = "ScanDeltaPtforSimTrackValid" + OutputPlotNameFix;
    OutputCanvas0->SaveAs(SaveName.c_str());

    TCanvas* OutputCanvas1 = new TCanvas("Canvas1", "Canvas1", 800, 600);
    OutputCanvas1->cd();

    //ScanEfficiencyforSimTrackValidGraph->SetTitle("ScanEfficienyforSimTrackValid");
    ScanEfficiencyforSimTrackValidGraph->GetXaxis()->SetTitle("simPt/Gev");
    ScanEfficiencyforSimTrackValidGraph->GetXaxis()->CenterTitle(1);
    ScanEfficiencyforSimTrackValidGraph->GetYaxis()->SetTitle("SeedingEfficiency %");
    ScanEfficiencyforSimTrackValidGraph->GetYaxis()->CenterTitle(1);
    ScanEfficiencyforSimTrackValidGraph->Draw("AP");
    string SaveName = "ScanEfficienyforSimTrackValid" + OutputPlotNameFix;
    OutputCanvas1->SaveAs(SaveName.c_str());

    //TCanvas* OutputCanvas2 = new TCanvas("Canvas2", "Canvas2", 800, 600);
    //OutputCanvas2->cd();


    //ScanChargeCheckforSimTrackValidGraph->SetTitle("ScanChargeCheckforSimTrackValid");
    ScanChargeCheckforSimTrackValidGraph->GetXaxis()->SetTitle("simPt/Gev");
    ScanChargeCheckforSimTrackValidGraph->GetXaxis()->CenterTitle(1);
    ScanChargeCheckforSimTrackValidGraph->GetYaxis()->SetTitle("InverseChargeRato %");
    ScanChargeCheckforSimTrackValidGraph->GetYaxis()->CenterTitle(1);
    ScanChargeCheckforSimTrackValidGraph->Draw("AP");
    string SaveName = "ScanChargeCheckforSimTrackValid" + OutputPlotNameFix;
    OutputCanvas1->SaveAs(SaveName.c_str());

    /*
       ScanChargeCheckforSimTrackValidHist->SetTitle("ScanChargeCheckforSimTrackValid");
       ScanChargeCheckforSimTrackValidHist->GetXaxis()->SetTitle("simPt/Gev");
       ScanChargeCheckforSimTrackValidHist->GetYaxis()->SetTitle("InverseChargeRato %");
       ScanChargeCheckforSimTrackValidHist->Draw("");
       string SaveName = "ScanChargeCheckforSimTrackValid" + OutputPlotNameFix;
       OutputCanvas0->SaveAs(SaveName.c_str());

       ScanEfficienyforSimTrackValidHist->SetTitle("ScanEfficienyforSimTrackValid");
       ScanEfficienyforSimTrackValidHist->GetXaxis()->SetTitle("simPt/Gev");
       ScanEfficienyforSimTrackValidHist->GetYaxis()->SetTitle("SeedingEfficiency %");
       ScanEfficienyforSimTrackValidHist->Draw("");
       string SaveName = "ScanEfficienyforSimTrackValid" + OutputPlotNameFix;
       OutputCanvas0->SaveAs(SaveName.c_str());
     */
}
