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

#define debug 0
#define recorde 1
#define F2PThresold 1.5
#define C2RThresold 20.0
#define PI 3.1415926
#define PhiBins 628
#define WorkPtRange 40.0
#define Chi2TH 5.0
#define FitOverRangeLimit 3
#define StatisticTH 100.0
#define isVertexConstraint 1
#define VertexLinkLimit 4

extern TStyle* gStyle;

using namespace std;

// old barrel layer index:0-5, new: 1-6
class BendingPhiIndexType {
    public:
        BendingPhiIndexType() {m[0] = 1; m[1] = 2; m[2] = 3; m[3] = 4;};
        BendingPhiIndexType(int a0, int a1, int a2, int a3) {m[0] = a0; m[1] = a1; m[2] = a2; m[3] = a3;};
        int m[4];
};

bool LessLayer(const unsigned int& a1, const unsigned int& a2) {
    bool isLess = false;
    if(a1 <= a2)
        isLess = true;

    return isLess;
}

class RPCBendingAnalyzer {
    public:
        RPCBendingAnalyzer();
        ~RPCBendingAnalyzer();
        void setParameters(string FileNamePara, string DrawOptionPara, int PatternTypePara, double PtScalePara, double BendingWiseCheckPara, double SegmentBendingCut0Para, double SegmentBendingCut1Para, double SegmentBendingCut2Para, double MaxBendingCutPara, bool applyFilterPara);
        void analyze(double thePhiC2R=-1.);
        void getEventRato(unsigned int IndexNumber, double theFitPtRange);
        void fitPtofPhi(unsigned int IndexNumber, string fitType, double thePhiC2R, double startPhiforMean, double endPhiforMean, double startPhiforSigma, double endPhiforSigma);
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
        int SampleLayer[10];
        double simHitBendingPhi[10][10];
        double simPtBendingPhi[10];
        double recHitBendingPhi[10][10];
        unsigned int RPCBendingLayer[10];
        int ClusterSize[10];
        int BX[10];
        double simPtatRef[10];

        //bool debug;
        TFile* File0;
        TTree* Tree0;
        string FileName;

        // PattentType: BarrelDouble->0 BarrelSingle1->1, EndcapSingle->2
        int PatternType;
        double PtScale;
        double SegmentBendingCut0;
        double SegmentBendingCut1;
        double SegmentBendingCut2;
        double MaxBendingCut;
        bool applyFilter;
        double BendingWiseCheck; // -1.: only reverse. 1.: only coverse. 0.: both

        double theRefPt;
        double simBendingPhiMax;
        double recBendingPhiMax;
        double recBendingPhiMin;
        double recBendingPhiMean;
        double BendingPhiVal0;
        double BendingPhiVal1;
        double BendingPhiVal2;
        int recBendingPhiMaxSign;
        int BendingPhiVal0Sign;
        int BendingPhiVal1Sign;
        int BendingPhiVal2Sign;
        double PhiC2R;
        double PhiF2P;
        double RatoThresoldC2R;
        double RatoThresoldF2P;
        std::vector<int> BendingPhiMaxFilter[15];
        //std::map<int, int> BendingPhiMaxFilterMap[15]; // map is much faster than vector in searching function

        string PatternFix;
        string CutFix;
        string BendingWiseFix;
        string FitFix;
        string FinalOutput;
        string OutputFix;

        unsigned int PatternIndex;
        std::vector<unsigned int> SampleLayerCollection;
        std::vector<BendingPhiIndexType> BendingPhiCollection;
        std::vector<double> BendingPhiResult;
        int validSample;
        TH2D* simHitBendingPhiHist[4][5][6][7];
        TH2D* recHitBendingPhiHist[4][5][6][7];
        TH2D* simHitBendingPhiMaxHist;
        TH2D* recHitBendingPhiMaxHist;
        TH2D* recHitBendingPhiMinHist;
        TH2D* recHitBendingPhiMeanHist;
        TH2D* recHitBendingPhiD0Hist;
        TH2D* recHitBendingPhiD1Hist;
        TH2D* recHitBendingPhiD2Hist;
        TH2D* recHitBendingPhiDSignHist;
        string theDrawOption;
};

RPCBendingAnalyzer::RPCBendingAnalyzer() {
}

RPCBendingAnalyzer::~RPCBendingAnalyzer() {
}

void RPCBendingAnalyzer::setParameters(string FileNamePara, string DrawOptionPara, int PatternTypePara, double PtScalePara, double BendingWiseCheckPara, double SegmentBendingCut0Para, double SegmentBendingCut1Para, double SegmentBendingCut2Para, double MaxBendingCutPara, bool applyFilterPara) {

    if(debug) cout << FileNamePara << endl;

    PatternType = PatternTypePara;
    PtScale = PtScalePara;
    SegmentBendingCut0 = SegmentBendingCut0Para;
    SegmentBendingCut1 = SegmentBendingCut1Para;
    SegmentBendingCut2 = SegmentBendingCut2Para;
    MaxBendingCut = MaxBendingCutPara;
    applyFilter = applyFilterPara;
    BendingWiseCheck = BendingWiseCheckPara;
    theDrawOption = DrawOptionPara;

    RatoThresoldF2P = F2PThresold;
    RatoThresoldC2R = C2RThresold;

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat("");

    if(debug) cout << "PatternType: " << PatternType << ", Draw option: " << theDrawOption << ", PtScale: " << PtScale << ", BendingWiseCheck: " << BendingWiseCheck << ", SegmentBendingCut0: " << SegmentBendingCut0 << ", SegmentBendingCut1: " << SegmentBendingCut1 << ", SegmentBendingCut2: " << SegmentBendingCut2 << ", MaxBendingCut: " << MaxBendingCut << ", applyFilter: " << applyFilter << endl;

    std::stringstream TempPtScale;
    TempPtScale << PtScale;
    string PtFix = TempPtScale.str() + "Gev_";

    if(PatternType%10==0)
        PatternFix = "oV_";
    if(PatternType%10==1)
        PatternFix = "wV_";

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

    CutFix = "";
    if(SegmentBendingCut0 == 0.)
        CutFix += "noCut0_";
    else {
        std::stringstream TempCut;
        TempCut << SegmentBendingCut0;
        string Cut = TempCut.str();
        CutFix += "Cut0>" + Cut + "_";
    }
    if(SegmentBendingCut1 == 0.)
        CutFix += "noCut1_";
    else {
        std::stringstream TempCut;
        TempCut << SegmentBendingCut1;
        string Cut = TempCut.str();
        CutFix += "Cut1>" + Cut + "_";
    }
    if(SegmentBendingCut2 == 0.)
        CutFix += "noCut2_";
    else {
        std::stringstream TempCut;
        TempCut << SegmentBendingCut2;
        string Cut = TempCut.str();
        CutFix += "Cut2>" + Cut + "_";
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
    Tree0->SetBranchAddress("SampleLayer", &SampleLayer);
    Tree0->SetBranchAddress("simHitBendingPhi", simHitBendingPhi);
    Tree0->SetBranchAddress("simPtBendingPhi", simPtBendingPhi);
    Tree0->SetBranchAddress("recHitBendingPhi", recHitBendingPhi);
    Tree0->SetBranchAddress("RPCBendingLayer", RPCBendingLayer);
    Tree0->SetBranchAddress("ClusterSize", ClusterSize);
    Tree0->SetBranchAddress("BX", BX);
    Tree0->SetBranchAddress("simPtatRef", simPtatRef);

    PatternIndex = -1;
    if((unsigned int)(PatternType/10) == 1234)
        PatternIndex = 0;
    if((unsigned int)(PatternType/10) == 1235)
        PatternIndex = 1;
    if((unsigned int)(PatternType/10) == 1236)
        PatternIndex = 2;
    if((unsigned int)(PatternType/10) == 1245)
        PatternIndex = 3;
    if((unsigned int)(PatternType/10) == 1246)
        PatternIndex = 4;
    if((unsigned int)(PatternType/10) == 1256)
        PatternIndex = 5;
    if((unsigned int)(PatternType/10) == 1345)
        PatternIndex = 6;
    if((unsigned int)(PatternType/10) == 1346)
        PatternIndex = 7;
    if((unsigned int)(PatternType/10) == 1356)
        PatternIndex = 8;
    if((unsigned int)(PatternType/10) == 1456)
        PatternIndex = 9;
    if((unsigned int)(PatternType/10) == 2345)
        PatternIndex = 10;
    if((unsigned int)(PatternType/10) == 2346)
        PatternIndex = 11;
    if((unsigned int)(PatternType/10) == 2356)
        PatternIndex = 12;
    if((unsigned int)(PatternType/10) == 2456)
        PatternIndex = 13;
    if((unsigned int)(PatternType/10) == 3456)
        PatternIndex = 14;
    if((unsigned int)(PatternType/10) == 135)
        PatternIndex = 15;
    if((unsigned int)(PatternType/10) == 136)
        PatternIndex = 16;
    if((unsigned int)(PatternType/10) == 146)
        PatternIndex = 17;
    if((unsigned int)(PatternType/10) == 246)
        PatternIndex = 18;
    if((unsigned int)(PatternType/10) == 156)
        PatternIndex = 19;
    if((unsigned int)(PatternType/10) == 256)
        PatternIndex = 20;    
    if((unsigned int)(PatternType/10) == 356)
        PatternIndex = 21;


    if(PatternIndex == -1) {
        if(debug) cout << "Pattern incorrect!" << endl;
        return;
    }

    unsigned int SampleLayer = (unsigned int)(PatternType/10);
    SampleLayerCollection.clear();
    while((unsigned int)(SampleLayer%10) > 0) {
        SampleLayerCollection.push_back((unsigned int)(SampleLayer%10));
        SampleLayer = (unsigned int)(SampleLayer/10);
    }
    sort(SampleLayerCollection.begin(), SampleLayerCollection.end(), LessLayer);


    if(SampleLayerCollection.size() < 3) {
        if(debug) cout << "not enough sample layers!" << endl;
        return;
    }

    BendingPhiCollection.clear();
    // barrel with vertex
    if(PatternType%10 == 1) {
        for(int i = 0; i < SampleLayerCollection.size()-1; i++)
            for(int j = i; j < SampleLayerCollection.size()-1; j++)
                for(int k = j; k < SampleLayerCollection.size()-1; k++)
                    for(int l = k+1; l < SampleLayerCollection.size(); l++) {
                        // we are in vertex constraint mode, so take vertex as much as possible
                        if(isVertexConstraint == 1)
                            if(i != j)
                                continue;
                        // only RB1/2 could link with vertex, since RB3 are close to RB4 and bring in large widthing
                        if(i == j && SampleLayerCollection[i] > VertexLinkLimit)
                            continue;
                        // close layers bring in large widthing
                        if((SampleLayerCollection[j] == (SampleLayerCollection[i]+1) && SampleLayerCollection[j] != 6) || (SampleLayerCollection[l] == (SampleLayerCollection[k]+1) && SampleLayerCollection[l] != 6))
                            continue;
                        BendingPhiCollection.push_back(BendingPhiIndexType(SampleLayerCollection[i], SampleLayerCollection[j], SampleLayerCollection[k], SampleLayerCollection[l]));
                    }
    }
    // barrel without vertex
    if(PatternType%10 == 0) {
        for(int i = 0; i < SampleLayerCollection.size()-2; i++)
            for(int j = i+1; j < SampleLayerCollection.size()-1; j++)
                for(int k = j; k < SampleLayerCollection.size()-1; k++)
                    for(int l = k+1; l < SampleLayerCollection.size(); l++) {
                        // close layers bring in large widthing
                        if(SampleLayerCollection[j] == (SampleLayerCollection[i]+1) || SampleLayerCollection[l] == (SampleLayerCollection[k]+1))
                            continue;
                        BendingPhiCollection.push_back(BendingPhiIndexType(SampleLayerCollection[i], SampleLayerCollection[j], SampleLayerCollection[k], SampleLayerCollection[l]));
                    }
    }

    // test
    if(1 == 2) {
        //int Filter0[] = {};
        //int Filter1[] = {3,4,5};

        //BendingPhiMaxFilter[0].assign(Filter0, Filter0+sizeof(Filter0)/sizeof(Filter0[0]));
        //BendingPhiMaxFilter[1].assign(Filter1, Filter1+sizeof(Filter1)/sizeof(Filter1[0]));

        //for(int i = 0; i < 8; i++)
        //for(int j = 0; j < BendingPhiMaxFilter[i].size(); j++)
        //BendingPhiMaxFilterMap[i][BendingPhiMaxFilter[i][j]] = -1;
    }

    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        //if(BendingPhiMaxFilterMap[PatternIndex][(int)Index] == -1)
        //continue;
        int i = BendingPhiCollection[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = BendingPhiCollection[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = BendingPhiCollection[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = BendingPhiCollection[Index].m[3];
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
    string recHitBendingPhiMeanHistName = FinalOutput + "recHitBendingPhiMeanHist";
    recHitBendingPhiMeanHist = new TH2D(recHitBendingPhiMeanHistName.c_str(), recHitBendingPhiMeanHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    string recHitBendingPhiD0HistName = FinalOutput + "recHitBendingPhiD0Hist";
    recHitBendingPhiD0Hist = new TH2D(recHitBendingPhiD0HistName.c_str(), recHitBendingPhiD0HistName.c_str(), PtScale*20, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    string recHitBendingPhiD1HistName = FinalOutput + "recHitBendingPhiD1Hist";
    recHitBendingPhiD1Hist = new TH2D(recHitBendingPhiD1HistName.c_str(), recHitBendingPhiD1HistName.c_str(), PtScale*20, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    string recHitBendingPhiD2HistName = FinalOutput + "recHitBendingPhiD2Hist";
    recHitBendingPhiD2Hist = new TH2D(recHitBendingPhiD2HistName.c_str(), recHitBendingPhiD2HistName.c_str(), PtScale*20, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    string recHitBendingPhiDSignHistName = FinalOutput + "recHitBendingPhiDSignHist";
    recHitBendingPhiDSignHist = new TH2D(recHitBendingPhiDSignHistName.c_str(), recHitBendingPhiDSignHistName.c_str(), PtScale*20, -1.*PtScale, PtScale, 5, -0.5, 4.5);
}

void RPCBendingAnalyzer::analyze(double thePhiC2R) {

    if(thePhiC2R > 0.)
        PhiC2R = thePhiC2R;

    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        //if(BendingPhiMaxFilterMap[PatternIndex][(int)Index] == -1)
        //continue;

        int i = BendingPhiCollection[Index].m[0];
        int j = BendingPhiCollection[Index].m[1];
        int k = BendingPhiCollection[Index].m[2];
        int l = BendingPhiCollection[Index].m[3];
        simHitBendingPhiHist[i][j][k][l]->Clear();
        recHitBendingPhiHist[i][j][k][l]->Clear();
    }
    simHitBendingPhiMaxHist->Clear();
    recHitBendingPhiMaxHist->Clear();
    recHitBendingPhiMinHist->Clear();
    recHitBendingPhiMeanHist->Clear();
    recHitBendingPhiDSignHist->Clear();

    int Nentries = Tree0->GetEntries();
    for(int i = 0; i < Nentries; i++) {
        Tree0->GetEntry(i);

        if(debug) cout << "PatternType: " << PatternType << ", simTrackMomentumPt: " << simTrackMomentumPt << endl;

        if(fabs(simTrackMomentumPt) > PtScale)
            continue;

        //if(fabs(simTrackMomentumEta) > 0.2)
        //continue;

        getBendingPhiMax();

        if(validSample == -1)
            continue;

        // Pattern Filter
        if(PatternType%10 == 1) {

            // choose different and far way RPCLayer for small widthing and less reverse bending in validating bendingPhi
            BendingPhiVal0 = getdPhi(recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[SampleLayerCollection.size()-2]], recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[0]]) * -1.;
            BendingPhiVal1 = getdPhi(recHitBendingPhi[SampleLayerCollection[1]][SampleLayerCollection[SampleLayerCollection.size()-1]], recHitBendingPhi[SampleLayerCollection[1]][SampleLayerCollection[1]]) * -1.;
            BendingPhiVal2 = getdPhi(recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[SampleLayerCollection.size()-1]], recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[0]]) * -1.;

            // filter small reverse bending if needed
            if(fabs(recBendingPhiMax) < MaxBendingCut || fabs(BendingPhiVal0) < SegmentBendingCut0 || fabs(BendingPhiVal1) < SegmentBendingCut1 || fabs(BendingPhiVal2) < SegmentBendingCut2) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            // for charge survey
            if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
                continue;

            recBendingPhiMaxSign = 1;
            BendingPhiVal0Sign = 1;
            BendingPhiVal1Sign = 1;
            BendingPhiVal2Sign = 1;
            if((double)simTrackCharge * recBendingPhiMax < 0.)
                recBendingPhiMaxSign = 0;
            if((double)simTrackCharge * BendingPhiVal0 < 0.)
                BendingPhiVal0Sign = 0;
            if((double)simTrackCharge * BendingPhiVal1 < 0.)
                BendingPhiVal1Sign = 0;
            if((double)simTrackCharge * BendingPhiVal2 < 0.)
                BendingPhiVal2Sign = 0;

            // check and correct reverse bending
            if(applyFilter == true) {
                int BendingPhiVal0FilterSign = (int)(BendingPhiVal0/fabs(BendingPhiVal0));
                int BendingPhiVal1FilterSign = (int)(BendingPhiVal1/fabs(BendingPhiVal1));
                int BendingPhiVal2FilterSign = (int)(BendingPhiVal2/fabs(BendingPhiVal2));
                int recBendingPhiMaxFilterSign = (int)(recBendingPhiMax/fabs(recBendingPhiMax));
                int BendingPhiSign = BendingPhiVal0FilterSign+BendingPhiVal1FilterSign+BendingPhiVal2FilterSign+recBendingPhiMaxFilterSign;
                if(abs(BendingPhiSign) < 2)
                    continue;
                /*
                   if((recBendingPhiMax * BendingPhiVal0 < 0. && SegmentBendingCut0 > 0.) || (recBendingPhiMax * BendingPhiVal1 < 0. && SegmentBendingCut1 > 0.) || (recBendingPhiMax * BendingPhiVal2 < 0. && SegmentBendingCut2 > 0.)) {

                   if(fabs(recBendingPhiMax) <= PhiC2R && SegmentBendingCut0 > 0. && SegmentBendingCut1 > 0. && recBendingPhiMax * BendingPhiVal0 < 0. && recBendingPhiMax * BendingPhiVal1 < 0.)
                   recBendingPhiMax *= -1.;
                   else {
                   if(debug) cout << "block by Filter." << endl;
                   continue;
                   }

                   int BendingPhiVal0FilterSign = 1;
                   int BendingPhiVal1FilterSign = 1;
                   int BendingPhiVal2FilterSign = 1;
                   if(recBendingPhiMax * BendingPhiVal0 < 0.)
                   BendingPhiVal0FilterSign = 0;
                   if(recBendingPhiMax * BendingPhiVal1 < 0.)
                   BendingPhiVal1FilterSign = 0;
                   if(recBendingPhiMax * BendingPhiVal2 < 0.)
                   BendingPhiVal2FilterSign = 0;
                   if((BendingPhiVal0FilterSign+BendingPhiVal1FilterSign+BendingPhiVal2FilterSign) <= 2)
                   continue;
                   }
                   */
            }
            theRefPt = simPtatRef[SampleLayerCollection[0]] * (double)simTrackCharge;
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
    recBendingPhiMean = 0.;
    double SampleCount = 0.;
    validSample = 1;

    double TempsimBendingPhi;
    double TemprecBendingPhi;
    BendingPhiResult.clear();
    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        //if(BendingPhiMaxFilterMap[PatternIndex][(int)Index] == -1)
        //continue;
        int i = BendingPhiCollection[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = BendingPhiCollection[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = BendingPhiCollection[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = BendingPhiCollection[Index].m[3];
        std::stringstream TempIndexL;
        TempIndexL << l;
        string IndexL = TempIndexL.str();

        if(ClusterSize[i] > 0 && ClusterSize[j] > 0 && ClusterSize[k] > 0 && ClusterSize[l] > 0) {
            TempsimBendingPhi = getdPhi(simHitBendingPhi[k][l], simHitBendingPhi[i][j]);
            TemprecBendingPhi = getdPhi(recHitBendingPhi[k][l], recHitBendingPhi[i][j]);
            if(i == j) {
                TempsimBendingPhi *= -1.;
                TemprecBendingPhi *= -1.;
            }
            BendingPhiResult.push_back(TemprecBendingPhi);
            if(fabs(simBendingPhiMax) < fabs(TempsimBendingPhi))
                simBendingPhiMax = TempsimBendingPhi;
            if(fabs(recBendingPhiMax) < fabs(TemprecBendingPhi))
                recBendingPhiMax = TemprecBendingPhi;
            if(fabs(recBendingPhiMin) > fabs(TemprecBendingPhi))
                recBendingPhiMin = TemprecBendingPhi;

            recBendingPhiMean += TemprecBendingPhi;
            SampleCount += 1.;
        }
        else
            validSample = -1;
    }
    recBendingPhiMean /= SampleCount;
}

void RPCBendingAnalyzer::fillBendingPhiHist() {

    if(debug) cout << "filling hists." << endl;

    if(validSample == -1)
        return;

    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        //if(BendingPhiMaxFilterMap[PatternIndex][(int)Index] == -1)
        //continue;
        int i = BendingPhiCollection[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = BendingPhiCollection[Index].m[1];
        std::stringstream TempIndexJ;
        TempIndexJ << j;
        string IndexJ = TempIndexJ.str();
        int k = BendingPhiCollection[Index].m[2];
        std::stringstream TempIndexK;
        TempIndexK << k;
        string IndexK = TempIndexK.str();
        int l = BendingPhiCollection[Index].m[3];
        std::stringstream TempIndexL;
        TempIndexL << l;
        string IndexL = TempIndexL.str();

        double simHitBendingPhiTemp = getdPhi(simHitBendingPhi[k][l], simHitBendingPhi[i][j]);
        double recHitBendingPhiTemp = getdPhi(recHitBendingPhi[k][l], recHitBendingPhi[i][j]);
        if(i == j) {
            simHitBendingPhiTemp *= -1.;
            recHitBendingPhiTemp *= -1.;
        }
        if(debug) cout << "simHitBendingPhiTemp: " << simHitBendingPhiTemp << ", recHitBendingPhiTemp: " << recHitBendingPhiTemp << endl;
        simHitBendingPhiHist[i][j][k][l]->Fill(theRefPt, simHitBendingPhiTemp);
        recHitBendingPhiHist[i][j][k][l]->Fill(theRefPt, recHitBendingPhiTemp);
    }
    simHitBendingPhiMaxHist->Fill(theRefPt, simBendingPhiMax);
    recHitBendingPhiMaxHist->Fill(theRefPt, recBendingPhiMax);
    recHitBendingPhiMinHist->Fill(theRefPt, recBendingPhiMin);
    recHitBendingPhiMeanHist->Fill(theRefPt, recBendingPhiMean);
    recHitBendingPhiDSignHist->Fill(theRefPt, BendingPhiVal0Sign+BendingPhiVal1Sign+BendingPhiVal2Sign+recBendingPhiMaxSign);

    if(recBendingPhiMax != 0.) {
        //cout << "test. recBendingPhiMax:" << recBendingPhiMax << ", BendingPhiVal0:" << BendingPhiVal0 << ", BendingPhiVal1:" << BendingPhiVal1 << ", BendingPhiVal2:" << BendingPhiVal2 << endl;
        recHitBendingPhiD0Hist->Fill(theRefPt, BendingPhiVal0/recBendingPhiMax);
        recHitBendingPhiD1Hist->Fill(theRefPt, BendingPhiVal1/recBendingPhiMax);
        recHitBendingPhiD2Hist->Fill(theRefPt, BendingPhiVal2/recBendingPhiMax);
    }
}

void RPCBendingAnalyzer::printHist() {

    TCanvas* BendingPhiCanvas = new TCanvas("BendingPhi", "BendingPhi", 800, 600);
    BendingPhiCanvas->cd();
    TPad* BendingPhiPad = new TPad("", "", 0, 0, 1, 1);
    BendingPhiPad->Draw();
    BendingPhiPad->cd();


    string SaveName;
    for(unsigned int Index = 0; Index < BendingPhiCollection.size(); Index++) {
        //if(BendingPhiMaxFilterMap[PatternIndex][(int)Index] == -1)
        //continue;
        int i = BendingPhiCollection[Index].m[0];
        std::stringstream TempIndexI;
        TempIndexI << i;
        string IndexI = TempIndexI.str();
        int j = BendingPhiCollection[Index].m[1];
        std::stringstream TempIndexJ;                   
        TempIndexJ << j;            
        string IndexJ = TempIndexJ.str();
        int k = BendingPhiCollection[Index].m[2];
        std::stringstream TempIndexK;                   
        TempIndexK << k;            
        string IndexK = TempIndexK.str();
        int l = BendingPhiCollection[Index].m[3];
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
    SaveName = FinalOutput + "recHitBendingMeanPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiMeanHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiMeanHist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiMeanHist->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiMeanHist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiMeanHist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());

    SaveName = FinalOutput + "recHitBendingD0Phi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiD0Hist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiD0Hist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiD0Hist->GetYaxis()->SetTitle("dPhi/Max");
    recHitBendingPhiD0Hist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiD0Hist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    SaveName = FinalOutput + "recHitBendingD1Phi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiD1Hist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiD1Hist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiD1Hist->GetYaxis()->SetTitle("dPhi/Max");
    recHitBendingPhiD1Hist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiD1Hist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    SaveName = FinalOutput + "recHitBendingD2Phi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiD2Hist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiD2Hist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiD2Hist->GetYaxis()->SetTitle("dPhi/Max");
    recHitBendingPhiD2Hist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiD2Hist->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    SaveName = FinalOutput + "recHitBendingDSignPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiDSignHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiDSignHist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiDSignHist->GetYaxis()->SetTitle("DSign");
    recHitBendingPhiDSignHist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiDSignHist->Draw("LEGO");
    BendingPhiCanvas->SaveAs(SaveName.c_str());

}

void RPCBendingAnalyzer::getEventRato(unsigned int IndexNumber, double theFitPtRange) {

    double FitPtRange = theFitPtRange;
    if(FitPtRange <= 0.)
        FitPtRange = WorkPtRange;

    unsigned int Index[4];
    for(unsigned int i = 0; i < 4; i++) {
        Index[i] = (unsigned int)(IndexNumber%10);
        IndexNumber = (unsigned int)(IndexNumber/10);
    }
    TH2D *SampleHist = recHitBendingPhiHist[Index[3]][Index[2]][Index[1]][Index[0]];
    
    string SampleHistName = SampleHist->GetName();

    string RatoC2RHistName = FinalOutput + SampleHistName + "_RatoC2R"; // Event of Coverse to Reverce bendingPhi
    TH1D* RatoC2RHist = new TH1D(RatoC2RHistName.c_str(), RatoC2RHistName.c_str(), 314, 0., PI/6.);
    string RatoF2PHistName = FinalOutput + SampleHistName + "_RatoF2P"; // Event of Fit region to Probe region
    TH1D* RatoF2PHist = new TH1D(RatoF2PHistName.c_str(), RatoF2PHistName.c_str(), 314, 0., PI/6.);

    int PtBinNumber = SampleHist->GetNbinsX();
    double PtLowEdge = SampleHist->GetXaxis()->GetBinLowEdge(1);
    double PtUpEdge = SampleHist->GetXaxis()->GetBinUpEdge(PtBinNumber);
    double PtBinWidth = SampleHist->GetXaxis()->GetBinWidth(1);
    int PhiBinNumber = SampleHist->GetNbinsY();
    double PhiLowEdge = SampleHist->GetYaxis()->GetBinLowEdge(1);
    double PhiUpEdge = SampleHist->GetYaxis()->GetBinUpEdge(PhiBinNumber);
    double PhiBinWidth = SampleHist->GetYaxis()->GetBinWidth(1);

    for(int PhiIndex = (int)(PhiBinNumber/2)+1; PhiIndex <= PhiBinNumber; PhiIndex++) {

        double PhiValue = (double)(PhiIndex) * PhiBinWidth + PhiLowEdge - PhiBinWidth * 0.5;
        if(PhiValue <= 0.)
            continue;

        std::stringstream TempIndex;
        TempIndex << PhiIndex;
        string PhiIndexName = TempIndex.str();

        string PtofPhiHistName = "PtofPhiHist_R" + PhiIndexName;
        TH1D* PtofPhiHist = SampleHist->ProjectionX(PtofPhiHistName.c_str(), PhiIndex, PhiIndex, "o");

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

void RPCBendingAnalyzer::fitPtofPhi(unsigned int IndexNumber, string fitType, double thePhiC2R, double startPhiforMean, double endPhiforMean, double startPhiforSigma, double endPhiforSigma) {

    unsigned int Index[4];
    for(unsigned int i = 0; i < 4; i++) {
        Index[i] = (unsigned int)(IndexNumber%10);
        IndexNumber = (unsigned int)(IndexNumber/10);
    }
    TH2D *SampleHist = recHitBendingPhiHist[Index[3]][Index[2]][Index[1]][Index[0]];

    gStyle->SetOptFit(1111);

    string SaveName;
    string SampleHistName = SampleHist->GetName();
    string FitHistName = FinalOutput + FitFix + SampleHistName;

    int PtBinNumber = SampleHist->GetNbinsX();
    double PtLowEdge = SampleHist->GetXaxis()->GetBinLowEdge(1);
    double PtUpEdge = SampleHist->GetXaxis()->GetBinUpEdge(PtBinNumber);
    double PtBinWidth = SampleHist->GetXaxis()->GetBinWidth(1);
    int PhiBinNumber = SampleHist->GetNbinsY();
    double PhiLowEdge = SampleHist->GetYaxis()->GetBinLowEdge(1);
    double PhiUpEdge = SampleHist->GetYaxis()->GetBinUpEdge(PhiBinNumber);
    double PhiBinWidth = SampleHist->GetYaxis()->GetBinWidth(1);
    if(debug) cout << "PtBinNumber: " << PtBinNumber << ", PtLowEdge: " << PtLowEdge << ", PtUpEdge: " << PtUpEdge << ", PtBinWidth: " << PtBinWidth << ". PhiBinNumber: " << PhiBinNumber << ", PhiLowEdge: " << PhiLowEdge << ", PhiUpEdge: " << PhiUpEdge << ", PhiBinWidth: " << PhiBinWidth << endl;

    string Phi2MeanPtHistName = FitHistName + "Phi2MeanPtHist";
    TH1D* Phi2MeanPtHist = new TH1D(Phi2MeanPtHistName.c_str(), Phi2MeanPtHistName.c_str(), PhiBins, 0., PI/6.);
    string Phi2RMSPtHistName = FitHistName + "Phi2RMSPtHist";
    TH1D* Phi2RMSPtHist = new TH1D(Phi2RMSPtHistName.c_str(), Phi2RMSPtHistName.c_str(), PhiBins, 0., PI/6.);
    string Phi2StatisticRatoHistName = FitHistName + "Phi2StatisticRatoHist";
    TH1D* Phi2StatisticRatoHist = new TH1D(Phi2StatisticRatoHistName.c_str(), Phi2StatisticRatoHistName.c_str(), 314, 0., PI/6.);

    double BendingPhiUpperLimit = PI/6.;
    bool PhiUpperLimitSet = false;
    double BendingPhiLowerLimit = PhiC2R; //PhiF2P;
    bool PhiLowerLimitSet = false;
    int outFit = 0;
    if(thePhiC2R >= 0.)
        BendingPhiLowerLimit = thePhiC2R;

    for(int index = (int)(PhiBinNumber/2)+1; index <= PhiBinNumber; index++) {
        double PhiValue = (double)(index) * PhiBinWidth + PhiLowEdge - PhiBinWidth * 0.5;

        // check only fit region
        if(PhiValue < BendingPhiLowerLimit)
            continue;

        std::stringstream TempIndex;
        TempIndex << index;
        string PhiIndex = TempIndex.str();


        string PtofPhiHistName = "PtofPhiHist_Test_" + PhiIndex;
        TH1D* PtofPhiHist = SampleHist->ProjectionX(PtofPhiHistName.c_str(), index, index, "o");
        double TotalEvent = 0.;
        double StartPtValue = 0.;
        double EndPtValue = WorkPtRange;//PtScale;
        bool isStartPtValueSet = false;
        bool isEndPtValueSet = false;
        int theHistMaxBin = PtofPhiHist->GetMaximumBin();
        double theHistMaxValue = PtofPhiHist->GetBinContent(theHistMaxBin);
        double theHistMaxBinCenter = PtofPhiHist->GetBinCenter(theHistMaxBin);
        for(int PtIndex = int((PtBinNumber/2)+1); PtIndex <= PtBinNumber; PtIndex++) {
            double PtValue = (double)(PtIndex) * PtBinWidth + PtLowEdge - PtBinWidth * 0.5;
            double TempValue = PtofPhiHist->GetBinContent(PtIndex);
            if((TempValue < (theHistMaxValue/50.0) || TempValue < 2) && PtIndex < theHistMaxBin && isStartPtValueSet == false) {
                StartPtValue = PtValue;
                //isStartPtValueSet = true;
            }
            if((TempValue < (theHistMaxValue/50.0) || TempValue < 2) && PtIndex > theHistMaxBin && isEndPtValueSet == false) {
                EndPtValue = PtValue;
                isEndPtValueSet = true;
                isStartPtValueSet = true;
            }
            TotalEvent += TempValue;
            if(debug) cout << "PtIndex: " << PtIndex << ", PtValue: " << PtValue << ", TempValue: " << TempValue << endl;
        }
        if(isStartPtValueSet == false)
            StartPtValue = 0.;        
        if(isEndPtValueSet == false)
            EndPtValue = PtScale;       

        // for this phi bin i we don't find event, so skip this one
        if(TotalEvent <= (double)StatisticTH)
            continue;

        string FitName = "PtofPhi_Fit_" + fitType;
        TF1* Fit0 = new TF1(FitName.c_str(), fitType.c_str(), StartPtValue, EndPtValue);
        TCanvas* PtofPhiCanvas = new TCanvas(PtofPhiHistName.c_str(), PtofPhiHistName.c_str(), 800, 600);
        PtofPhiCanvas->cd();
        TPad* PtofPhiPad = new TPad(PtofPhiHistName.c_str(), PtofPhiHistName.c_str(), 0, 0, 1, 1);
        PtofPhiPad->Draw();
        PtofPhiPad->cd();
        /*
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
        */
        int NewBinNumber = 5; // rebin to 5x
        string FinalPtofPhiHistName = "PtofPhiHist_Final_" + PhiIndex;
        TH1D* FinalPtofPhiHist = (TH1D*)PtofPhiHist->Rebin(NewBinNumber, FinalPtofPhiHistName.c_str());
        int NewPtBinNumber = FinalPtofPhiHist->GetNbinsX();
        double NewPtLowEdge = FinalPtofPhiHist->GetBinLowEdge(1);
        double NewPtBinWidth = FinalPtofPhiHist->GetBinWidth(1);
        /*
        bool NewLargeStatistic = false;
        double EventNumberInsideRegion = 0.;
        double EventNumberOutsideRegion = 0.;
        TotalEvent = 0;
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
        
        if(NewLargeStatistic == false) {
           EndPtValue /= 2.;
           Fit0->SetRange(StartPtValue, EndPtValue);
        }

        */
        double EventC = 0.;
        double EventR = 0.;
        isStartPtValueSet = false;
        isEndPtValueSet = false;
        int theFinalHistMaxBin = FinalPtofPhiHist->GetMaximumBin();
        double theFinalHistMaxValue = FinalPtofPhiHist->GetBinContent(theFinalHistMaxBin);
        double theFinalHistMaxBinCenter = FinalPtofPhiHist->GetBinCenter(theFinalHistMaxBin);
        for(int PtIndex = 1; PtIndex <= NewPtBinNumber; PtIndex++) {
            double PtValue = (double)(PtIndex) * NewPtBinWidth + NewPtLowEdge - NewPtBinWidth * 0.5;
            double TempValue = FinalPtofPhiHist->GetBinContent(PtIndex);

            if(PtValue < 0.)
                EventR += TempValue;
            else
                EventC += TempValue;

            if((TempValue < (theFinalHistMaxValue/50.0) || TempValue < 2) && PtIndex < theFinalHistMaxBin && isStartPtValueSet == false) {
                StartPtValue = PtValue;
                //isStartPtValueSet = true;
            }
            if((TempValue < (theFinalHistMaxValue/50.0) || TempValue < 2) && PtIndex > theFinalHistMaxBin && isEndPtValueSet == false) {
                EndPtValue = PtValue;
                isEndPtValueSet = true;
                isStartPtValueSet = true;
            }
            //TotalEvent += TempValue;
            if(debug) cout << "PtIndex: " << PtIndex << ", PtValue: " << PtValue << ", TempValue: " << TempValue << endl;
        }
        double Rato = EventC / (EventR+0.01);
        if(isStartPtValueSet == false)
            StartPtValue = 0.;
        if(isEndPtValueSet == false)
            EndPtValue = PtScale;

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

        if(FinalMean < StartPtValue || FinalMean >= EndPtValue || fabs(FinalMean-theFinalHistMaxBinCenter) > 10. ) { //|| FinalReducedChi3 > Chi2TH) { 
            FinalMean = theFinalHistMaxBinCenter;
            FinalSigma = FinalMean / 5.;
            isFinalCorrected = true;
            if(debug) cout << "Correct Final: " << FinalMean << ", " << FinalSigma << endl;
        }

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
        if(isFinalCorrected == false) {
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
        Phi2StatisticRatoHist->GetYaxis()->SetTitle("Rato");
        Phi2StatisticRatoHist->GetYaxis()->CenterTitle(1);
        Phi2StatisticRatoHist->SetAxisRange(0., 20., "Y");
        Phi2StatisticRatoHist->Draw("");
        SaveName = FitHistName + "Phi2StatisticRato" + OutputFix;
        Phi2StatisticRatoCanvas->SaveAs(SaveName.c_str());

    }
