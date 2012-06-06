#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include <algorithm>

#define debug 1
#define recorde 1
#define F2PThresold 1.5
#define C2RThresold 4.0
#define PI 3.1415926
#define PhiBins 628
#define WorkPtRange 40.0
#define Chi2TH 5.0
#define FitOverRangeLimit 3
#define StatisticTH 100.0

extern TStyle* gStyle;

using namespace std;

class BendingPhiIndexType {
    public:
        BendingPhiIndexType() {m[0] = 0; m[1] = 0; m[2] = 0; m[3] = 0;};
        BendingPhiIndexType(int a0, int a1, int a2, int a3) {m[0] = a0; m[1] = a1; m[2] = a2; m[3] = a3;};
        int m[4];
};

class RPCBendingAnalyzer {
    public:
        RPCBendingAnalyzer();
        ~RPCBendingAnalyzer();
        void setParameters(string FileNamePara, string DrawOptionPara, int PatternTypePara, double PtScalePara, double BendingWiseCheckPara, double SegmentBendingCutPara, double MaxBendingCutPara, bool applyFilterPara);
        void analyze(double thePhiC2R=-1.);
        void getEventRato(double theFitPtRange);
        void fitPtofPhi(string fitType, double thePhiF2P, double startPhiforMean, double endPhiforMean, double startPhiforSigma, double endPhiforSigma);
    private:
        void getBendingPhiMax();
        double getdPhi(double Phi1, double Phi0);
        void fillBendingPhiHist();
        void printHist();

        unsigned int EventNumber;
        int simTrackId;
        int simTrackType;
        double simTrackMomentumPt;
        double simTrackMomentumEta;
        int simTrackCharge;
        int simTrackvalid;
        double simHitBendingPhi[6][6];
        double simPtBendingPhi[6];
        double recHitBendingPhi[6][6];
        unsigned int RPCBendingLayer[6];
        int ClusterSize[6];
        int BX[6];
        double simPtatRef[6];

        //bool debug;
        TFile* File0;
        TTree* Tree0;
        string FileName;

        // PattentType: BarrelDouble->0 BarrelSingle1->1, EndcapSingle->2
        int PatternType;
        double PtScale;
        double SegmentBendingCut;
        double MaxBendingCut;
        bool applyFilter;
        double BendingWiseCheck; // -1.: only reverse. 1.: only coverse. 0.: both

        double theRefPt;
        double simBendingPhiMax;
        double recBendingPhiMax;
        double recBendingPhiMin;
        double PhiC2R;
        double PhiF2P;
        double RatoThresoldC2R;
        double RatoThresoldF2P;
        std::vector<int> BendingPhiMaxFilter[8];
        std::map<int, int> BendingPhiMaxFilterMap[8]; // map is much faster than vector in searching function

        string PatternFix;
        string CutFix;
        string BendingWiseFix;
        string FitFix;
        string FinalOutput;
        string OutputFix;

        std::vector<BendingPhiIndexType> HistFilter;
        TH2D* simHitBendingPhiHist[3][4][4][6];
        TH2D* recHitBendingPhiHist[3][4][4][6];
        TH2D* simHitBendingPhiMaxHist;
        TH2D* recHitBendingPhiMaxHist;
        TH2D* recHitBendingPhiMinHist;
        string theDrawOption;
};

RPCBendingAnalyzer::RPCBendingAnalyzer() {
}

RPCBendingAnalyzer::~RPCBendingAnalyzer() {
}

void RPCBendingAnalyzer::setParameters(string FileNamePara, string DrawOptionPara, int PatternTypePara, double PtScalePara, double BendingWiseCheckPara, double SegmentBendingCutPara, double MaxBendingCutPara, bool applyFilterPara) {

    if(debug) cout << FileNamePara << endl;

    PatternType = PatternTypePara;
    PtScale = PtScalePara;
    SegmentBendingCut = SegmentBendingCutPara;
    MaxBendingCut = MaxBendingCutPara;
    applyFilter = applyFilterPara;
    BendingWiseCheck = BendingWiseCheckPara;
    theDrawOption = DrawOptionPara;

    RatoThresoldF2P = F2PThresold;
    RatoThresoldC2R = C2RThresold;

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat("");

    if(debug) cout << "PatternType: " << PatternType << ", Draw option: " << theDrawOption << ", PtScale: " << PtScale << ", BendingWiseCheck: " << BendingWiseCheck << ", SegmentBendingCut: " << SegmentBendingCut << ", MaxBendingCut: " << MaxBendingCut << ", applyFilter: " << applyFilter << endl;

    std::stringstream TempPtScale;
    TempPtScale << PtScale;
    string PtFix = TempPtScale.str() + "Gev_";

    if(PatternType%10==0)
        PatternFix = "BD_";
    if(PatternType%10==1)
        PatternFix = "BS_";
    if(PatternType%10==2)
        PatternFix = "ES_";

    PatternFix += "Algo";
    std::stringstream AlgorithmChoice;
    AlgorithmChoice << (int)(PatternType/10);
    PatternFix += AlgorithmChoice.str();
    PatternFix += "_";

    if(BendingWiseCheck < 0.)
        BendingWiseFix = "Reverse_";
    if(BendingWiseCheck > 0.)
        BendingWiseFix = "Coverse_";
    if(BendingWiseCheck == 0.)
        BendingWiseFix = "Both_";

    if(SegmentBendingCut == 0.)
        CutFix = "noCut1234_";
    else {
        std::stringstream TempCut;
        TempCut << SegmentBendingCut;
        string Cut = TempCut.str();
        CutFix = "Cut1234>" + Cut + "_";
    }

    if(MaxBendingCut == 0.)
        CutFix += "noCutMax_";
    else {
        std::stringstream TempCut;
        TempCut << MaxBendingCut;
        string Cut = TempCut.str();
        CutFix += "CutMax>" + Cut + "_";
    }

    if(applyFilter == false)
        CutFix += "noFilter_";
    else
        CutFix += "Filter_";

    FinalOutput = theDrawOption + "_" + PtFix + PatternFix + BendingWiseFix + CutFix;
    OutputFix = ".eps";

    //string theFileName = "rfio:/castor/cern.ch/user/h/hyteng/rpcseed/validation/Pt3.0-PtScale.0Gev_Eta-1.0-1.0_recBending362/recBA_key_merge/" + FileName;
    FileName = FileNamePara;
    string theFileName = FileName;

    File0 = (TFile*)TFile::Open(theFileName.c_str());
    Tree0 = (TTree*)File0->Get("ExTree");

    Tree0->SetBranchAddress("EventNumber", &EventNumber);
    Tree0->SetBranchAddress("simTrackId", &simTrackId);
    Tree0->SetBranchAddress("simTrackType", &simTrackType);
    Tree0->SetBranchAddress("simTrackMomentumPt", &simTrackMomentumPt);
    Tree0->SetBranchAddress("simTrackMomentumEta", &simTrackMomentumEta);
    Tree0->SetBranchAddress("simTrackCharge", &simTrackCharge);
    Tree0->SetBranchAddress("simTrackvalid", &simTrackvalid);
    Tree0->SetBranchAddress("simHitBendingPhi", simHitBendingPhi);
    Tree0->SetBranchAddress("simPtBendingPhi", simPtBendingPhi);
    Tree0->SetBranchAddress("recHitBendingPhi", recHitBendingPhi);
    Tree0->SetBranchAddress("RPCBendingLayer", RPCBendingLayer);
    Tree0->SetBranchAddress("ClusterSize", ClusterSize);
    Tree0->SetBranchAddress("BX", BX);
    Tree0->SetBranchAddress("simPtatRef", simPtatRef);

    HistFilter.clear();
    if(PatternType%10 == 0) {
        HistFilter.push_back(BendingPhiIndexType(0,1,2,3));
        HistFilter.push_back(BendingPhiIndexType(0,1,1,4));
        HistFilter.push_back(BendingPhiIndexType(0,1,1,5));
        HistFilter.push_back(BendingPhiIndexType(2,3,3,4));
        HistFilter.push_back(BendingPhiIndexType(2,3,3,5));
        HistFilter.push_back(BendingPhiIndexType(0,0,0,1));
        HistFilter.push_back(BendingPhiIndexType(2,2,2,3));

        int Filter0[] = {0,2,4};
        int Filter1[] = {2,4,5,6}; // kick out 0
        int Filter2[] = {0,1,3};
        int Filter3[] = {1,3,5,6}; // kick out 0
        int Filter4[] = {0,1,2,3,4,5,6}; // fullset

        BendingPhiMaxFilter[0].assign(Filter0, Filter0+sizeof(Filter0)/sizeof(Filter0[0]));
        BendingPhiMaxFilter[1].assign(Filter1, Filter1+sizeof(Filter1)/sizeof(Filter1[0]));
        BendingPhiMaxFilter[2].assign(Filter2, Filter2+sizeof(Filter2)/sizeof(Filter2[0]));
        BendingPhiMaxFilter[3].assign(Filter3, Filter3+sizeof(Filter3)/sizeof(Filter3[0]));
        BendingPhiMaxFilter[4].assign(Filter4, Filter4+sizeof(Filter4)/sizeof(Filter4[0]));

        for(int i = 0; i < 8; i++)
            for(int j = 0; j < BendingPhiMaxFilter[i].size(); j++)
                BendingPhiMaxFilterMap[i][BendingPhiMaxFilter[i][j]] = -1;
    }
    if(PatternType%10 == 1) {
        HistFilter.push_back(BendingPhiIndexType(0,1,1,2));
        HistFilter.push_back(BendingPhiIndexType(0,1,1,3));
        HistFilter.push_back(BendingPhiIndexType(0,1,1,4));
        HistFilter.push_back(BendingPhiIndexType(0,1,1,5));
        HistFilter.push_back(BendingPhiIndexType(2,3,3,4));
        HistFilter.push_back(BendingPhiIndexType(2,3,3,5));
        HistFilter.push_back(BendingPhiIndexType(0,0,0,1));
        HistFilter.push_back(BendingPhiIndexType(2,2,2,3));

        int Filter0[] = {2,3};
        int Filter1[] = {2,3,6};
        int Filter2[] = {0,1,3};
        int Filter3[] = {0,1,3,6};
        int Filter4[] = {0,1,2};
        int Filter5[] = {0,1,2,6};
        int Filter6[] = {4,5};
        int Filter7[] = {4,5,7};

        BendingPhiMaxFilter[0].assign(Filter0, Filter0+sizeof(Filter0)/sizeof(Filter0[0]));
        BendingPhiMaxFilter[1].assign(Filter1, Filter1+sizeof(Filter1)/sizeof(Filter1[0]));
        BendingPhiMaxFilter[2].assign(Filter2, Filter2+sizeof(Filter2)/sizeof(Filter2[0]));
        BendingPhiMaxFilter[3].assign(Filter3, Filter3+sizeof(Filter3)/sizeof(Filter3[0]));
        BendingPhiMaxFilter[4].assign(Filter4, Filter4+sizeof(Filter4)/sizeof(Filter4[0]));
        BendingPhiMaxFilter[5].assign(Filter5, Filter5+sizeof(Filter5)/sizeof(Filter5[0]));
        BendingPhiMaxFilter[6].assign(Filter6, Filter6+sizeof(Filter6)/sizeof(Filter6[0]));
        BendingPhiMaxFilter[7].assign(Filter7, Filter7+sizeof(Filter7)/sizeof(Filter7[0]));

        for(int i = 0; i < 8; i++) 
            for(int j = 0; j < BendingPhiMaxFilter[i].size(); j++) 
                BendingPhiMaxFilterMap[i][BendingPhiMaxFilter[i][j]] = -1;

    }
    if(PatternType%10 == 2) {
        HistFilter.push_back(BendingPhiIndexType(0,3,3,4));
        HistFilter.push_back(BendingPhiIndexType(0,3,3,5));
        HistFilter.push_back(BendingPhiIndexType(0,4,4,5));
        HistFilter.push_back(BendingPhiIndexType(0,0,0,3));
        HistFilter.push_back(BendingPhiIndexType(0,0,0,4));
        HistFilter.push_back(BendingPhiIndexType(0,0,0,5));

        int Filter0[] = {0,1,2};
        int Filter1[] = {3,4,5};

        BendingPhiMaxFilter[0].assign(Filter0, Filter0+sizeof(Filter0)/sizeof(Filter0[0]));
        BendingPhiMaxFilter[1].assign(Filter1, Filter1+sizeof(Filter1)/sizeof(Filter1[0]));

        for(int i = 0; i < 8; i++)
            for(int j = 0; j < BendingPhiMaxFilter[i].size(); j++)
                BendingPhiMaxFilterMap[i][BendingPhiMaxFilter[i][j]] = -1;
    }

    for(unsigned int Index = 0; Index < HistFilter.size(); Index++) {
        if(BendingPhiMaxFilterMap[PatternType/10][(int)Index] != -1)
            continue;
        int i = HistFilter[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = HistFilter[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = HistFilter[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = HistFilter[Index].m[3];
        std::stringstream TempIndexL;
        TempIndexL << l;
        string IndexL = TempIndexL.str();

        string simHitBendingPhiHistName = FinalOutput + "simHitBendingPhiHist" + IndexI + IndexJ + IndexK + IndexL;
        simHitBendingPhiHist[i][j][k][l] = new TH2D(simHitBendingPhiHistName.c_str(), simHitBendingPhiHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins*10, -1.*PI/6., PI/6.);
        string recHitBendingPhiHistName = FinalOutput + "recHitBendingPhiHist" + IndexI + IndexJ + IndexK + IndexL;
        recHitBendingPhiHist[i][j][k][l] = new TH2D(recHitBendingPhiHistName.c_str(), recHitBendingPhiHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    }

    string simHitBendingPhiMaxHistName = FinalOutput + "simHitBendingPhiMaxHist";
    simHitBendingPhiMaxHist = new TH2D(simHitBendingPhiMaxHistName.c_str(), simHitBendingPhiMaxHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins*10, -1.*PI/6., PI/6.);
    string recHitBendingPhiMaxHistName = FinalOutput + "recHitBendingPhiMaxHist";
    recHitBendingPhiMaxHist = new TH2D(recHitBendingPhiMaxHistName.c_str(), recHitBendingPhiMaxHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    string recHitBendingPhiMinHistName = FinalOutput + "recHitBendingPhiMinHist";
    recHitBendingPhiMinHist = new TH2D(recHitBendingPhiMinHistName.c_str(), recHitBendingPhiMinHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
}

void RPCBendingAnalyzer::analyze(double thePhiC2R) {

    if(thePhiC2R > 0.)
        PhiC2R = thePhiC2R;

    for(unsigned int Index = 0; Index < HistFilter.size(); Index++) {
        if(BendingPhiMaxFilterMap[PatternType/10][(int)Index] != -1)
            continue;

        int i = HistFilter[Index].m[0];
        int j = HistFilter[Index].m[1];
        int k = HistFilter[Index].m[2];
        int l = HistFilter[Index].m[3];
        simHitBendingPhiHist[i][j][k][l]->Clear();
        recHitBendingPhiHist[i][j][k][l]->Clear();
    }
    simHitBendingPhiMaxHist->Clear();
    recHitBendingPhiMaxHist->Clear();
    recHitBendingPhiMinHist->Clear();

    int Nentries = Tree0->GetEntries();
    for(int i = 0; i < Nentries; i++) {
        Tree0->GetEntry(i);

        if(debug) cout << "PatternType: " << PatternType << ", simTrackMomentumPt: " << simTrackMomentumPt << endl;

        if(fabs(simTrackMomentumPt) > PtScale)
            continue;

        //if(fabs(simTrackMomentumEta) > 0.2)
        //continue;

        getBendingPhiMax();

        // Pattern Filter
        if(PatternType == 40) {
            //if(fabs(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1])) < SegmentBendingCut)
            //continue;

            //if(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1])*simTrackCharge > 0.)
            //continue;

            //if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. || recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0. || recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1]) < 0.)
            //continue;

            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 0 || PatternType == 10) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            //if(recBendingPhiMax*(double)simTrackCharge > 0.)
            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) <= PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. || recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0.) {
                        if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. && recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0.)// || (recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. && recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1]) > 0.))
                            recBendingPhiMax *= -1.; // correct the charge
                        else {
                            if(debug) cout << "block by Filter." << endl;
                            continue;
                        }
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 20 || PatternType == 30) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[0][1])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) <= PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. || recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0. ) {
                        if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0. && recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0.)
                            recBendingPhiMax *= -1.; // correct the charge
                        else {
                            if(debug) cout << "block by Filter." << endl;
                            continue;
                        }
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 1 || PatternType == 11) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) < PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0.) {
                        if(debug) cout << "block by Filter." << endl;
                        continue;
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 21 || PatternType == 31) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) < PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0.) {
                        if(debug) cout << "block by Filter." << endl;
                        continue;
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 41 || PatternType == 51) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) < PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[0][1], recHitBendingPhi[0][0]) > 0.) {
                        if(debug) cout << "block by Filter." << endl;
                        continue;
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        if(PatternType == 61 || PatternType == 71) {

            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) < PhiC2R)
                    if(recBendingPhiMax * getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2]) > 0.) {
                        if(debug) cout << "block by Filter." << endl;
                        continue;
                    }
            theRefPt = simPtatRef[2] * (double)simTrackCharge;
        }

        if(PatternType == 2 || PatternType == 12) {
            if(fabs(recBendingPhiMax) < MaxBendingCut && fabs(getdPhi(recHitBendingPhi[2][3], recHitBendingPhi[2][2])) < SegmentBendingCut) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            if(applyFilter == true)
                if(fabs(recBendingPhiMax) < PhiC2R)
                    if((recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][3]) > 0.) || (recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][4]) > 0.) || (recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][5]) > 0.)) {
                        if((recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][3]) > 0.) && (recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][4]) > 0.) && (recBendingPhiMax * getdPhi(recHitBendingPhi[0][0], recHitBendingPhi[0][5]) > 0.))
                            recBendingPhiMax *= -1.; // correct the charge
                        else {
                            if(debug) cout << "block by Filter." << endl;
                            continue;
                        }
                    }
            theRefPt = simPtatRef[0] * (double)simTrackCharge;
        }

        fillBendingPhiHist();
    }

    if(debug) cout << "done filling." << endl;
    printHist();
}

double RPCBendingAnalyzer::getdPhi(double Phi1, double Phi0) {
    double dPhi = Phi1 - Phi0;
    if(fabs(dPhi) > PI)
        dPhi = 2. * PI - dPhi;
    if(debug) cout << "Phi1: " << Phi1 << ", Phi0: " << Phi0 << ", dPhi: " << dPhi << endl;
    return dPhi;
}

void RPCBendingAnalyzer::getBendingPhiMax() {

    simBendingPhiMax = 0.;
    recBendingPhiMax = 0.;
    recBendingPhiMin = 10.;

    double TempsimBendingPhi;
    double TemprecBendingPhi;
    for(unsigned int Index = 0; Index < HistFilter.size(); Index++) {
        if(BendingPhiMaxFilterMap[PatternType/10][(int)Index] != -1)
            continue;
        int i = HistFilter[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = HistFilter[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = HistFilter[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = HistFilter[Index].m[3];
        std::stringstream TempIndexL;
        TempIndexL << l;
        string IndexL = TempIndexL.str();

        TempsimBendingPhi = getdPhi(simHitBendingPhi[k][l], simHitBendingPhi[i][j]);
        TemprecBendingPhi = getdPhi(recHitBendingPhi[k][l], recHitBendingPhi[i][j]);
        if(i == j) {
            TempsimBendingPhi *= -1.;
            TemprecBendingPhi *= -1.;
        }
        if(fabs(simBendingPhiMax) < fabs(TempsimBendingPhi))
            simBendingPhiMax = TempsimBendingPhi;
        if(fabs(recBendingPhiMax) < fabs(TemprecBendingPhi))
            recBendingPhiMax = TemprecBendingPhi;
        if(fabs(recBendingPhiMin) > fabs(TemprecBendingPhi))
            recBendingPhiMin = TemprecBendingPhi;
    }
}

void RPCBendingAnalyzer::fillBendingPhiHist() {

    if(debug) cout << "filling hists." << endl;

    for(unsigned int Index = 0; Index < HistFilter.size(); Index++) {
        if(BendingPhiMaxFilterMap[PatternType/10][(int)Index] != -1)
            continue;
        int i = HistFilter[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = HistFilter[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = HistFilter[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = HistFilter[Index].m[3];
        std::stringstream TempIndexL;
        TempIndexL << l;
        string IndexL = TempIndexL.str();

        double simHitBendingPhiTemp = getdPhi(simHitBendingPhi[k][l], simHitBendingPhi[i][j]);
        double recHitBendingPhiTemp = getdPhi(recHitBendingPhi[k][l], recHitBendingPhi[i][j]);
        if(debug) cout << "simHitBendingPhiTemp: " << simHitBendingPhiTemp << ", recHitBendingPhiTemp: " << recHitBendingPhiTemp << endl;
        simHitBendingPhiHist[i][j][k][l]->Fill(theRefPt, simHitBendingPhiTemp);
        recHitBendingPhiHist[i][j][k][l]->Fill(theRefPt, recHitBendingPhiTemp);
    }
    simHitBendingPhiMaxHist->Fill(theRefPt, simBendingPhiMax);
    recHitBendingPhiMaxHist->Fill(theRefPt, recBendingPhiMax);
    recHitBendingPhiMinHist->Fill(theRefPt, recBendingPhiMin);
}

void RPCBendingAnalyzer::printHist() {

    TCanvas* BendingPhiCanvas = new TCanvas("BendingPhi", "BendingPhi", 800, 600);
    BendingPhiCanvas->cd();
    TPad* BendingPhiPad = new TPad("", "", 0, 0, 1, 1);
    BendingPhiPad->Draw();
    BendingPhiPad->cd();


    string SaveName;
    for(unsigned int Index = 0; Index < HistFilter.size(); Index++) {
        if(BendingPhiMaxFilterMap[PatternType/10][(int)Index] != -1)
            //if(find(BendingPhiMaxFilter[PatternType/10].begin(), BendingPhiMaxFilter[PatternType/10].end(), (int)Index) == BendingPhiMaxFilter[PatternType/10].end())
            continue;
        int i = HistFilter[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = HistFilter[Index].m[1];
        std::stringstream TempIndexJ;                   
        TempIndexJ << j;            
        string IndexJ = TempIndexJ.str();
        int k = HistFilter[Index].m[2];
        std::stringstream TempIndexK;                   
        TempIndexK << k;            
        string IndexK = TempIndexK.str();
        int l = HistFilter[Index].m[3];
        std::stringstream TempIndexL;                   
        TempIndexL << l;            
        string IndexL = TempIndexL.str();
        if(debug) cout << "BendingPhiIndex: " << i << j << k << l << ". simIntegral: " << simHitBendingPhiHist[i][j][k][l]->Integral() << ", recIntegral: " << recHitBendingPhiHist[i][j][k][l]->Integral() << endl;
        SaveName = FinalOutput + "simHitBendingPhi" + IndexI + IndexJ + IndexK + IndexL + OutputFix;
        if(debug) cout << SaveName << endl;
        simHitBendingPhiHist[i][j][k][l]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
        simHitBendingPhiHist[i][j][k][l]->GetXaxis()->CenterTitle(1);
        simHitBendingPhiHist[i][j][k][l]->GetYaxis()->SetTitle("Bending in Phi");
        simHitBendingPhiHist[i][j][k][l]->GetYaxis()->CenterTitle(1);
        simHitBendingPhiHist[i][j][k][l]->SetStats(0);
        simHitBendingPhiHist[i][j][k][l]->Draw(theDrawOption.c_str());
        BendingPhiCanvas->SaveAs(SaveName.c_str());
        SaveName = FinalOutput + "recHitBendingPhi" + IndexI + IndexJ + IndexK + IndexL + OutputFix;
        if(debug) cout << SaveName << endl;
        recHitBendingPhiHist[i][j][k][l]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
        recHitBendingPhiHist[i][j][k][l]->GetXaxis()->CenterTitle(1);
        recHitBendingPhiHist[i][j][k][l]->GetYaxis()->SetTitle("Bending in Phi");
        recHitBendingPhiHist[i][j][k][l]->GetYaxis()->CenterTitle(1);
        recHitBendingPhiHist[i][j][k][l]->Draw(theDrawOption.c_str());
        BendingPhiCanvas->SaveAs(SaveName.c_str());
    }

    std::stringstream TempIndex;
    TempIndex << PatternType;
    string Index = TempIndex.str();
    SaveName = FinalOutput + "simHitBendingMaxPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    simHitBendingPhiMaxHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    simHitBendingPhiMaxHist->GetXaxis()->CenterTitle(1);
    simHitBendingPhiMaxHist->GetYaxis()->SetTitle("Bending in Phi");
    simHitBendingPhiMaxHist->GetYaxis()->CenterTitle(1);
    simHitBendingPhiMaxHist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    SaveName = FinalOutput + "recHitBendingMaxPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiMaxHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiMaxHist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiMaxHist->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiMaxHist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiMaxHist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    SaveName = FinalOutput + "recHitBendingMinPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiMinHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiMinHist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiMinHist->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiMinHist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiMinHist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
}

void RPCBendingAnalyzer::getEventRato(double theFitPtRange) {

    double FitPtRange = theFitPtRange;
    if(FitPtRange < 0.)
        FitPtRange = WorkPtRange;

    string RatoC2RHistName = FinalOutput + "RatoC2R"; // Event of Coverse to Reverce bendingPhi
    TH1D* RatoC2RHist = new TH1D(RatoC2RHistName.c_str(), RatoC2RHistName.c_str(), 314, 0., PI/6.);
    string RatoF2PHistName = FinalOutput + "RatoF2P"; // Event of Fit region to Probe region
    TH1D* RatoF2PHist = new TH1D(RatoF2PHistName.c_str(), RatoF2PHistName.c_str(), 314, 0., PI/6.);

    int PtBinNumber = recHitBendingPhiMaxHist->GetNbinsX();
    double PtLowEdge = recHitBendingPhiMaxHist->GetXaxis()->GetBinLowEdge(1);
    double PtUpEdge = recHitBendingPhiMaxHist->GetXaxis()->GetBinUpEdge(PtBinNumber);
    double PtBinWidth = recHitBendingPhiMaxHist->GetXaxis()->GetBinWidth(1);
    int PhiBinNumber = recHitBendingPhiMaxHist->GetNbinsY();
    double PhiLowEdge = recHitBendingPhiMaxHist->GetYaxis()->GetBinLowEdge(1);
    double PhiUpEdge = recHitBendingPhiMaxHist->GetYaxis()->GetBinUpEdge(PhiBinNumber);
    double PhiBinWidth = recHitBendingPhiMaxHist->GetYaxis()->GetBinWidth(1);

    for(int PhiIndex = (int)(PhiBinNumber/2)+1; PhiIndex <= PhiBinNumber; PhiIndex++) {

        double PhiValue = (double)(PhiIndex) * PhiBinWidth + PhiLowEdge - PhiBinWidth * 0.5;
        if(PhiValue <= 0.)
            continue;

        std::stringstream TempIndex;
        TempIndex << PhiIndex;
        string PhiIndexName = TempIndex.str();

        string PtofPhiHistName = "PtofPhiHist_R" + PhiIndexName;
        TH1D* PtofPhiHist = recHitBendingPhiMaxHist->ProjectionX(PtofPhiHistName.c_str(), PhiIndex, PhiIndex, "o");

        int PtBinNumber = PtofPhiHist->GetNbinsX();
        double PtLowEdge = PtofPhiHist->GetBinLowEdge(1);
        double PtBinWidth = PtofPhiHist->GetBinWidth(1);
        double EventNumberInsideFitRegion = 0.;
        double EventNumberOutsideFitRegion = 0.;
        double EventNumberCoverse = 0.;
        double EventNumberReverse = 0.;

        for(int PtIndex = 1; PtIndex <= PtBinNumber; PtIndex++) {

            double PtValue = (double)(PtIndex) * PtBinWidth + PtLowEdge - PtBinWidth * 0.5;
            double TempEventNumber = PtofPhiHist->GetBinContent(PtIndex);

            if(PtValue > FitPtRange || PtValue < 0.)
                EventNumberOutsideFitRegion += TempEventNumber;
            else
                EventNumberInsideFitRegion += TempEventNumber;

            if(PtValue > 0.)
                EventNumberCoverse += TempEventNumber;
            else
                EventNumberReverse += TempEventNumber;
        }
        delete PtofPhiHist;
        EventNumberReverse += 0.1;
        RatoC2RHist->SetBinContent((PhiIndex-(int)(PhiBinNumber/2)), EventNumberCoverse/EventNumberReverse);
        if(EventNumberCoverse/EventNumberReverse <= RatoThresoldC2R && PhiValue <= 0.12)
            PhiC2R = PhiValue;
        EventNumberOutsideFitRegion += 0.1;
        RatoF2PHist->SetBinContent((PhiIndex-(int)(PhiBinNumber/2)), EventNumberInsideFitRegion/EventNumberOutsideFitRegion);
        if(EventNumberInsideFitRegion/EventNumberOutsideFitRegion <= RatoThresoldF2P && PhiValue <= 0.12)
            PhiF2P = PhiValue;

        if(recorde) cout << "PhiC2R: " << PhiC2R << ", PhiF2P: " << PhiF2P << endl;
    }

    string SaveName;
    TCanvas* Rato2PhiCanvas = new TCanvas("Rato2Phi", "Rato2Phi", 800, 600);
    Rato2PhiCanvas->cd();
    TPad* Rato2PhiPad = new TPad("Rato2Phi", "Rato2Phi", 0, 0, 1, 1);
    Rato2PhiPad->Draw();
    Rato2PhiPad->cd();
    RatoC2RHist->GetXaxis()->SetTitle("max bending in Phi");
    RatoC2RHist->GetXaxis()->CenterTitle(1);
    RatoC2RHist->GetYaxis()->SetTitle("rato of C2R");
    RatoC2RHist->GetYaxis()->CenterTitle(1);
    RatoC2RHist->SetAxisRange(0., 10., "Y");
    RatoC2RHist->Draw();
    SaveName = RatoC2RHistName + ".png";
    Rato2PhiCanvas->SaveAs(SaveName.c_str());
    RatoF2PHist->GetXaxis()->SetTitle("max bending in Phi");
    RatoF2PHist->GetXaxis()->CenterTitle(1);
    RatoF2PHist->GetYaxis()->SetTitle("rato of F2P");
    RatoF2PHist->GetYaxis()->CenterTitle(1);
    RatoF2PHist->SetAxisRange(0., 10., "Y");
    RatoF2PHist->Draw();
    SaveName = RatoF2PHistName + ".png";
    Rato2PhiCanvas->SaveAs(SaveName.c_str());

    delete RatoC2RHist;
    delete RatoF2PHist;
}

void RPCBendingAnalyzer::fitPtofPhi(string fitType, double thePhiF2P, double startPhiforMean, double endPhiforMean, double startPhiforSigma, double endPhiforSigma) {

    string SaveName;
    string FitHistName = FinalOutput + FitFix + "recHitBendingPhiMax_";

    int PtBinNumber = recHitBendingPhiMaxHist->GetNbinsX();
    double PtLowEdge = recHitBendingPhiMaxHist->GetXaxis()->GetBinLowEdge(1);
    double PtUpEdge = recHitBendingPhiMaxHist->GetXaxis()->GetBinUpEdge(PtBinNumber);
    double PtBinWidth = recHitBendingPhiMaxHist->GetXaxis()->GetBinWidth(1);
    int PhiBinNumber = recHitBendingPhiMaxHist->GetNbinsY();
    double PhiLowEdge = recHitBendingPhiMaxHist->GetYaxis()->GetBinLowEdge(1);
    double PhiUpEdge = recHitBendingPhiMaxHist->GetYaxis()->GetBinUpEdge(PhiBinNumber);
    double PhiBinWidth = recHitBendingPhiMaxHist->GetYaxis()->GetBinWidth(1);
    if(debug) cout << "PtBinNumber: " << PtBinNumber << ", PtLowEdge: " << PtLowEdge << ", PtUpEdge: " << PtUpEdge << ", PtBinWidth: " << PtBinWidth << ". PhiBinNumber: " << PhiBinNumber << ", PhiLowEdge: " << PhiLowEdge << ", PhiUpEdge: " << PhiUpEdge << ", PhiBinWidth: " << PhiBinWidth << endl;

    string Phi2MeanPtHistName = FitHistName + "Phi2MeanPtHist";
    TH1D* Phi2MeanPtHist = new TH1D(Phi2MeanPtHistName.c_str(), Phi2MeanPtHistName.c_str(), PhiBins, 0., PI/6.);
    string Phi2RMSPtHistName = FitHistName + "Phi2RMSPtHist";
    TH1D* Phi2RMSPtHist = new TH1D(Phi2RMSPtHistName.c_str(), Phi2RMSPtHistName.c_str(), PhiBins, 0., PI/6.);
    string Phi2StatisticRatoHistName = FitHistName + "Phi2StatisticRatoHist";
    TH1D* Phi2StatisticRatoHist = new TH1D(Phi2StatisticRatoHistName.c_str(), Phi2StatisticRatoHistName.c_str(), 314, 0., PI/6.);

    double BendingPhiUpperLimit = PI/6.;
    bool PhiUpperLimitSet = false;
    double BendingPhiLowerLimit = PhiF2P;
    bool PhiLowerLimitSet = false;
    int outFit = 0;
    if(thePhiF2P > 0.)
        BendingPhiLowerLimit = thePhiF2P;

    for(int index = (int)(PhiBinNumber/2)+1; index <= PhiBinNumber; index++) {
        double PhiValue = (double)(index) * PhiBinWidth + PhiLowEdge - PhiBinWidth * 0.5;

        if(PhiValue < BendingPhiLowerLimit)
            continue;

        std::stringstream TempIndex;
        TempIndex << index;
        string PhiIndex = TempIndex.str();


        string PtofPhiHistName = "PtofPhiHist_Test_" + PhiIndex;
        TH1D* PtofPhiHist = recHitBendingPhiMaxHist->ProjectionX(PtofPhiHistName.c_str(), index, index, "o");
        double TotalEvent = 0.;
        double StartPtValue = 0.;
        double EndPtValue = WorkPtRange;
        for(int PtIndex = (int)(PtBinNumber/2)+1; PtIndex <= PtBinNumber; PtIndex++) {

            double PtValue = (double)(PtIndex) * PtBinWidth + PtLowEdge - PtBinWidth * 0.5;
            double TempValue = PtofPhiHist->GetBinContent(PtIndex);
            TotalEvent += TempValue;
            if(debug) cout << "PtIndex: " << PtIndex << ", PtValue: " << PtValue << ", TempValue: " << TempValue << endl;
        }

        // for this phi bin i we don't find event, so skip this one
        if(TotalEvent <= (double)StatisticTH)
            continue;

        string FitName = "PtofPhi_Fit_" + fitType;
        TF1* Fit0 = new TF1(FitName.c_str(), fitType.c_str(), StartPtValue, EndPtValue);
        double theHistMaxValue = PtofPhiHist->GetBinContent(PtofPhiHist->GetMaximumBin());
        double theHistMaxBinCenter = PtofPhiHist->GetBinCenter(PtofPhiHist->GetMaximumBin());
        Fit0->SetParameter(0, theHistMaxValue);
        Fit0->SetParameter(1, theHistMaxBinCenter);
        PtofPhiHist->Fit(Fit0, "N", "", StartPtValue, EndPtValue);

        gStyle->SetOptFit(0);
        TCanvas* PtofPhiCanvas = new TCanvas(PtofPhiHistName.c_str(), PtofPhiHistName.c_str(), 800, 600);
        PtofPhiCanvas->cd();
        TPad* PtofPhiPad = new TPad(PtofPhiHistName.c_str(), PtofPhiHistName.c_str(), 0, 0, 1, 1);
        PtofPhiPad->Draw();
        PtofPhiPad->cd();
        PtofPhiHist->Draw("");
        Fit0->SetLineStyle(2);
        Fit0->SetLineColor(4);
        Fit0->Draw("same");
        SaveName = FitHistName + PtofPhiHistName + OutputFix;
        PtofPhiCanvas->SaveAs(SaveName.c_str());

        int NewBinNumber = 5;
        string FinalPtofPhiHistName = "PtofPhiHist_Final_" + PhiIndex;
        TH1D* FinalPtofPhiHist = (TH1D*)PtofPhiHist->Rebin(NewBinNumber, FinalPtofPhiHistName.c_str());
        int NewPtBinNumber = FinalPtofPhiHist->GetNbinsX();
        double NewPtLowEdge = FinalPtofPhiHist->GetBinLowEdge(1);
        double NewPtBinWidth = FinalPtofPhiHist->GetBinWidth(1);
        bool NewLargeStatistic = false;
        double EventNumberInsideRegion = 0.;
        double EventNumberOutsideRegion = 0.;

        for(int PtIndex = 1; PtIndex <= NewPtBinNumber; PtIndex++) {

            double PtValue = (double)(PtIndex) * NewPtBinWidth + NewPtLowEdge - NewPtBinWidth * 0.5;
            double TempValue = FinalPtofPhiHist->GetBinContent(PtIndex);

            if(PtValue > WorkPtRange || PtValue < StartPtValue)
                EventNumberOutsideRegion += TempValue;
            else
                EventNumberInsideRegion += TempValue;

            if(debug) cout << "PtIndex: " << PtIndex << ", PtValue: " << PtValue << ", TempValue: " << TempValue << endl;
        }
        if(EventNumberOutsideRegion == 0.)
            EventNumberOutsideRegion += 0.1;
        double Rato = EventNumberInsideRegion/EventNumberOutsideRegion;
        if(debug) cout << "NewLargeStatistic: " << NewLargeStatistic << ", EventNumberInsideRegion: " << EventNumberInsideRegion << ", EventNumberOutsideRegion: " << EventNumberOutsideRegion << ", Rato: " << Rato << endl;
        std::stringstream RegionStatisticRatoTemp;
        RegionStatisticRatoTemp << Rato ;
        string RegionStatisticRato = RegionStatisticRatoTemp.str();
        string RegionStatisticInfo = "RegionStatisticInfo: " + RegionStatisticRato;
        /*
           if(NewLargeStatistic == false) {
           EndPtValue /= 2.;
           Fit0->SetRange(StartPtValue, EndPtValue);
           }
           */
        double theFinalHistMaxValue = FinalPtofPhiHist->GetBinContent(FinalPtofPhiHist->GetMaximumBin());
        double theFinalHistMaxBinCenter = FinalPtofPhiHist->GetBinCenter(FinalPtofPhiHist->GetMaximumBin());
        Fit0->SetParameter(0, theFinalHistMaxValue);
        Fit0->SetParameter(1, theFinalHistMaxBinCenter);
        FinalPtofPhiHist->Fit(Fit0, "", "", StartPtValue, EndPtValue);

        double FinalChiSquare = Fit0->GetChisquare();
        double FinalMean = Fit0->GetParameter(1);
        double FinalSigma = Fit0->GetParameter(2);
        double FinalReducedChi2 = Fit0->GetChisquare() / Fit0->GetNDF();
        if(debug) cout << "FinalChiSquare: " << FinalChiSquare << ", FinalMean: " << FinalMean << ", FinalSigma:" << FinalSigma << endl;
        // when fitting fail Mean will over PtScale range, while Sigma keep in same level. we confirm this from Fit plots
        bool isFinalCorrected = false;
        if(FinalMean < 0. || FinalMean > WorkPtRange ) { //|| FinalReducedChi3 > Chi2TH) { 
            FinalMean = theFinalHistMaxBinCenter;
            FinalSigma = FinalMean / 5.;
            isFinalCorrected = true;
            if(debug) cout << "Correct Final: " << FinalMean << ", " << FinalSigma << endl;
        }

        gStyle->SetOptFit(1111);
        PtofPhiCanvas->cd();
        PtofPhiPad->cd();
        FinalPtofPhiHist->Draw("");
        //PtofPhiPad->Update();
        //TPaveStats *St = (TPaveStats*)FinalPtofPhiHist->FindObject("stats");
        //St->InsertText(RegionStatisticInfo.c_str());
        //St->Paint("");
        //PtofPhiPad->Update();
        //FinalPtofPhiHist->SetStats(1);
        Fit0->SetLineStyle(2);
        Fit0->SetLineColor(3);
        Fit0->Draw("same");
        SaveName = FitHistName + FinalPtofPhiHistName + OutputFix;
        PtofPhiCanvas->SaveAs(SaveName.c_str());

        delete PtofPhiHist;
        delete FinalPtofPhiHist;
        delete Fit0;
        delete PtofPhiPad;
        delete PtofPhiCanvas;

        double MeanPt = FinalMean;
        double SigmaPt = FinalSigma;
        if(debug) cout << "index: " << index << ", PhiValue: " << PhiValue << ", MeanPt: " << MeanPt << ", SigmaPt: " << SigmaPt << endl;
        int thePhiIndex = Phi2MeanPtHist->FindBin(PhiValue);
        if(Rato >= 1.5 && isFinalCorrected == false) {
            outFit = 0;
            Phi2MeanPtHist->SetBinContent(thePhiIndex, MeanPt);
            Phi2RMSPtHist->SetBinContent(thePhiIndex, SigmaPt);
            if(PhiValue < 0.1 && PhiLowerLimitSet == false) {
                BendingPhiLowerLimit = PhiValue;
                PhiLowerLimitSet = true;
            }
            if(PhiValue > 0.1 && PhiUpperLimitSet == false)
                BendingPhiUpperLimit = PhiValue;
        }
        else {
            outFit++;
            if(PhiValue < 0.1 && outFit >= FitOverRangeLimit)
                PhiLowerLimitSet = false;
            if(PhiUpperLimitSet == false && outFit >= FitOverRangeLimit)
                PhiUpperLimitSet = true;
        }
        Phi2StatisticRatoHist->SetBinContent(thePhiIndex, Rato);
        }
        if(recorde) std::cout << "BendingPhiLowerLimit: " << BendingPhiLowerLimit << ", BendingPhiUpperLimit: " << BendingPhiUpperLimit << endl;

        if(startPhiforMean <= 0.)
            startPhiforMean = BendingPhiLowerLimit;
        if(endPhiforMean <= 0.)
            endPhiforMean = BendingPhiUpperLimit;
        if(startPhiforSigma <= 0.)
            startPhiforSigma = BendingPhiLowerLimit;
        if(endPhiforSigma <= 0.)
            endPhiforSigma = BendingPhiUpperLimit;
        if(recorde) std::cout << "startPhiforMean: " << startPhiforMean << ", endPhiforMean: " << endPhiforMean << ", startPhiforSigma: " << startPhiforSigma << ", endPhiforSigma: " << endPhiforSigma << endl;

        gStyle->SetOptFit(0111);
        double para[3];
        TF1 *PtValueFunction = new TF1("fitPtValue", "[0]*x*x+[1]*x+[2]", startPhiforMean, endPhiforMean);
        //PtValueFunction->SetParameters(a, b, c);
        Phi2MeanPtHist->Fit(PtValueFunction, "R");
        PtValueFunction->GetParameters(para);
        if(recorde) cout << "Fitted para: " << para[0] << ", " << para[1] << ", " << para[2] << endl;

        TF1 *PtErrorFunction = new TF1("fitPtError", "[0]*x*x+[1]*x+[2]", startPhiforSigma, endPhiforSigma);
        //PtErrorFunction->SetParameters(a, b, c);
        Phi2RMSPtHist->Fit(PtErrorFunction, "R");
        PtErrorFunction->GetParameters(para);
        if(recorde) cout << "Fitted para: " << para[0] << ", " << para[1] << ", " << para[2] << endl;


        TCanvas* Phi2MeanPtCanvas = new TCanvas("Phi2MeanPt", "Phi2MeanPt", 800, 600);
        Phi2MeanPtCanvas->cd();
        TPad* Phi2MeanPtPad = new TPad("Phi2MeanPt", "Phi2MeanPt", 0, 0, 1, 1);
        Phi2MeanPtPad->Draw();
        Phi2MeanPtPad->cd();
        Phi2MeanPtHist->GetXaxis()->SetTitle("max bending in Phi");
        Phi2MeanPtHist->GetXaxis()->CenterTitle(1);
        Phi2MeanPtHist->GetYaxis()->SetTitle("mean Pt estimation Gev");
        Phi2MeanPtHist->GetYaxis()->CenterTitle(1);
        Phi2MeanPtHist->Draw("");
        PtValueFunction->SetLineStyle(2);
        PtValueFunction->SetLineColor(3);
        PtValueFunction->Draw("same");
        SaveName = FitHistName + "Phi2MeanPt" + OutputFix;
        Phi2MeanPtCanvas->SaveAs(SaveName.c_str());

        TCanvas* Phi2RMSPtCanvas = new TCanvas("Phi2RMSPt", "Phi2RMSPt", 800, 600);
        Phi2RMSPtCanvas->cd();
        TPad* Phi2RMSPtPad = new TPad("Phi2RMSPt", "Phi2RMSPt", 0, 0, 1, 1);
        Phi2RMSPtPad->Draw();
        Phi2RMSPtPad->cd();
        Phi2RMSPtHist->GetXaxis()->SetTitle("max bending in Phi");
        Phi2RMSPtHist->GetXaxis()->CenterTitle(1);
        Phi2RMSPtHist->GetYaxis()->SetTitle("sigma Pt estimation Gev");
        Phi2RMSPtHist->GetYaxis()->CenterTitle(1);
        Phi2RMSPtHist->Draw("");
        PtErrorFunction->SetLineStyle(2);
        PtErrorFunction->SetLineColor(3);
        PtErrorFunction->Draw("same");
        SaveName = FitHistName + "Phi2RMSPt" + OutputFix;
        Phi2RMSPtCanvas->SaveAs(SaveName.c_str());

        TCanvas* Phi2StatisticRatoCanvas = new TCanvas("Phi2StatisticRato", "Phi2StatisticRato", 800, 600);
        Phi2StatisticRatoCanvas->cd();
        TPad* Phi2StatisticRatoPad = new TPad("Phi2StatisticRato", "Phi2StatisticRato", 0, 0, 1, 1);
        Phi2StatisticRatoPad->Draw();
        Phi2StatisticRatoPad->cd();
        Phi2StatisticRatoHist->GetXaxis()->SetTitle("max bending in Phi");
        Phi2StatisticRatoHist->GetXaxis()->CenterTitle(1);
        Phi2StatisticRatoHist->GetYaxis()->SetTitle("mean Pt estimation Gev");
        Phi2StatisticRatoHist->GetYaxis()->CenterTitle(1);
        Phi2StatisticRatoHist->Draw("");
        SaveName = FitHistName + "Phi2StatisticRato" + OutputFix;
        Phi2StatisticRatoCanvas->SaveAs(SaveName.c_str());
    }
