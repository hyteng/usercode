#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include <algorithm>

#define debug 1
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
//#define VertexLinkLimit 5

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
        void runSample(double ms0, double me0, double ss0, double se0, double ms1, double me1, double ss1, double se1, double ms2, double me2, double ss2, double se2, double ms3, double me3, double ss3, double se3);

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
        double VertexLinkLimit;
        double theRefPt;
        double simBendingPhiMax;
        double recBendingPhiMax;
        double recBendingPhiMin;
        double recBendingPhiMean;
        double BendingPhiVal[4];
        int BendingPhiValSign[4];
        int SimBendingPhiValFilterSign[4];
        int BendingPhiValFilterSign[4];
        int SimBendingPhiFilterSign;
        int BendingPhiFilterSign;
        double PhiC2R;
        double PhiF2P;
        double RatoThresoldC2R;
        double RatoThresoldF2P;
        //std::vector<int> BendingPhiMaxFilter[15];
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
        TH2F* simHitBendingPhiHist[6][6][6][7];
        TH2F* recHitBendingPhiHist[6][6][6][7];
        TH2F* simHitBendingPhiMaxHist;
        TH2F* recHitBendingPhiMaxHist;
        TH2F* recHitBendingPhiMinHist;
        TH2F* recHitBendingPhiMeanHist;
        TH2F* recHitBendingPhiDHist[3]; // D0,D1,D2 w.r.t Max
        TH2F* recHitBendingPhiSHist[4]; // D3+Max1
        TH2F* recHitBendingPhiBHist[2];
        TH2F* recHitBendingPhiDSignHist;
        string theDrawOption;


        double ParaMean[4][3];
        double ParaSigma[4][3];
        double StartPhiMean[4];
        double EndPhiMean[4];
        double StartPhiSigma[4];
        double EndPhiSigma[4];
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
    OutputFix = ".png";

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

    if(SampleLayerCollection.size() > 3)
        VertexLinkLimit = 4;
    else
        VertexLinkLimit = 5;


    if(SampleLayerCollection.size() < 3) {
        if(debug) cout << "not enough sample layers!" << endl;
        return;
    }

    BendingPhiCollection.clear();
    // barrel with vertex
    //if(PatternType%10 == 1) {
        for(int i = 0; i < SampleLayerCollection.size()-1; i++)
            for(int j = i; j < SampleLayerCollection.size()-1; j++)
                for(int k = j; k < SampleLayerCollection.size()-1; k++)
                    for(int l = k+1; l < SampleLayerCollection.size(); l++) {
                        // we are in vertex constraint mode, so take vertex as much as possible
                        if(isVertexConstraint == 1)
                            if(i != j || j != k || (SampleLayerCollection.size() == 4 && l == k+1))
                                continue;
                        // only RB1/2 could link with vertex, since RB3 are close to RB4 and bring in large widthing
                        if(i == j && SampleLayerCollection[i] > VertexLinkLimit)
                            continue;
                        // close layers bring in large widthing
                        if((SampleLayerCollection[j] == (SampleLayerCollection[i]+1) && SampleLayerCollection[j] != 6) || (SampleLayerCollection[l] == (SampleLayerCollection[k]+1) && SampleLayerCollection[l] != 6))
                            continue;

                        if(debug) cout << "sampling BendingIndex with Vertex: " << SampleLayerCollection[i]  << SampleLayerCollection[j] << SampleLayerCollection[k] << SampleLayerCollection[l] << endl;
                        BendingPhiCollection.push_back(BendingPhiIndexType(SampleLayerCollection[i], SampleLayerCollection[j], SampleLayerCollection[k], SampleLayerCollection[l]));
                        //delete simHitBendingPhiHist[SampleLayerCollection[i]][SampleLayerCollection[j]][SampleLayerCollection[k]][SampleLayerCollection[l]];
                        //delete recHitBendingPhiHist[SampleLayerCollection[i]][SampleLayerCollection[j]][SampleLayerCollection[k]][SampleLayerCollection[l]];
                    }
    //}
    /*
    // barrel without vertex
    if(PatternType%10 == 0) {
        for(int i = 0; i < SampleLayerCollection.size()-2; i++)
            for(int j = i+1; j < SampleLayerCollection.size()-1; j++)
                for(int k = j; k < SampleLayerCollection.size()-1; k++)
                    for(int l = k+1; l < SampleLayerCollection.size(); l++) {
                        // close layers bring in large widthing
                        if(SampleLayerCollection[j] == (SampleLayerCollection[i]+1) || SampleLayerCollection[l] == (SampleLayerCollection[k]+1))
                            continue;

                        if(debug) cout << "sampling BendingIndex without Vertex: " << SampleLayerCollection[i]  << SampleLayerCollection[j] << SampleLayerCollection[k] << SampleLayerCollection[l] << endl;
                        BendingPhiCollection.push_back(BendingPhiIndexType(SampleLayerCollection[i], SampleLayerCollection[j], SampleLayerCollection[k], SampleLayerCollection[l]));
                    }
    }
    */
    /*
    // test
    for(int i = 0; i < 5; i++)
            for(int j = 0; j < 5; j++)
                for(int k = 0; k < 5; k++)
                    for(int l = 0; l < 6; l++) {
                        delete simHitBendingPhiHist[i][j][k][l];
                        delete recHitBendingPhiHist[i][j][k][l];
                    }
    */
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

        delete simHitBendingPhiHist[i][j][k][l];
        delete recHitBendingPhiHist[i][j][k][l];

        if(debug) cout << "allocating histgrame with BendingIndex: " << i << j << k << l << endl;
        string simHitBendingPhiHistName = FinalOutput + "simHitBendingPhiHist" + IndexI + IndexJ + IndexK + IndexL;
        simHitBendingPhiHist[i][j][k][l] = new TH2F(simHitBendingPhiHistName.c_str(), simHitBendingPhiHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
        string recHitBendingPhiHistName = FinalOutput + "recHitBendingPhiHist" + IndexI + IndexJ + IndexK + IndexL;
        recHitBendingPhiHist[i][j][k][l] = new TH2F(recHitBendingPhiHistName.c_str(), recHitBendingPhiHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    }

    string simHitBendingPhiMaxHistName = FinalOutput + "simHitBendingPhiMaxHist";
    simHitBendingPhiMaxHist = new TH2F(simHitBendingPhiMaxHistName.c_str(), simHitBendingPhiMaxHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    string recHitBendingPhiMaxHistName = FinalOutput + "recHitBendingPhiMaxHist";
    recHitBendingPhiMaxHist = new TH2F(recHitBendingPhiMaxHistName.c_str(), recHitBendingPhiMaxHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    string recHitBendingPhiMinHistName = FinalOutput + "recHitBendingPhiMinHist";
    recHitBendingPhiMinHist = new TH2F(recHitBendingPhiMinHistName.c_str(), recHitBendingPhiMinHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    string recHitBendingPhiMeanHistName = FinalOutput + "recHitBendingPhiMeanHist";
    recHitBendingPhiMeanHist = new TH2F(recHitBendingPhiMeanHistName.c_str(), recHitBendingPhiMeanHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);

    for(int i = 0; i < 3; i++) {
        std::stringstream TempIndexI;
        TempIndexI << i;
        string recHitBendingPhiDHistName = FinalOutput + "recHitBendingPhiD" + TempIndexI.str() + "Hist";
        recHitBendingPhiDHist[i] = new TH2F(recHitBendingPhiDHistName.c_str(), recHitBendingPhiDHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 300, -1.5, 1.5);
        string recHitBendingPhiSHistName = FinalOutput + "recHitBendingPhiS" + TempIndexI.str() + "Hist";
        recHitBendingPhiSHist[i] = new TH2F(recHitBendingPhiSHistName.c_str(), recHitBendingPhiSHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    }
    string recHitBendingPhiSHistName = FinalOutput + "recHitBendingPhiS3Hist";
    recHitBendingPhiSHist[3] = new TH2F(recHitBendingPhiSHistName.c_str(), recHitBendingPhiSHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, PhiBins, -1.*PI/6., PI/6.);
    //recHitBendingPhiSHist[3] = new TH2F(recHitBendingPhiSHistName.c_str(), recHitBendingPhiSHistName.c_str(), PtScale*2, -1., 1., PhiBins, -1.*PI/6., PI/6.);
    /*
    string recHitBendingPhiD0HistName = FinalOutput + "recHitBendingPhiD0Hist";
    recHitBendingPhiD0Hist = new TH2F(recHitBendingPhiD0HistName.c_str(), recHitBendingPhiD0HistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    string recHitBendingPhiD1HistName = FinalOutput + "recHitBendingPhiD1Hist";
    recHitBendingPhiD1Hist = new TH2F(recHitBendingPhiD1HistName.c_str(), recHitBendingPhiD1HistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    string recHitBendingPhiD2HistName = FinalOutput + "recHitBendingPhiD2Hist";
    recHitBendingPhiD2Hist = new TH2F(recHitBendingPhiD2HistName.c_str(), recHitBendingPhiD2HistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 500, -2.5, 2.5);
    */
    string recHitBendingPhiDSignHistName = FinalOutput + "recHitBendingPhiDSignHist";
    recHitBendingPhiDSignHist = new TH2F(recHitBendingPhiDSignHistName.c_str(), recHitBendingPhiDSignHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 5, -0.5, 4.5);

    string recHitBendingPhiBHistName = FinalOutput + "recHitBendingPhiB10"  + "Hist";
    recHitBendingPhiBHist[0] = new TH2F(recHitBendingPhiBHistName.c_str(), recHitBendingPhiBHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 300, -3., 3.);
    recHitBendingPhiBHistName = FinalOutput + "recHitBendingPhiB20"  + "Hist";
    recHitBendingPhiBHist[1] = new TH2F(recHitBendingPhiBHistName.c_str(), recHitBendingPhiBHistName.c_str(), PtScale*10, -1.*PtScale, PtScale, 300, -3., 3.);
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

        if(fabs(simTrackMomentumEta) > 1.0)
            continue;

        getBendingPhiMax();

        if(validSample == -1)
            continue;

        // Pattern Filter
        if(PatternType%10 == 1) {

            // choose different and far way RPCLayer for small widthing and less reverse bending in validating bendingPhi
            BendingPhiVal[0] = getdPhi(recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[SampleLayerCollection.size()-2]], recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[0]]) * -1.;
            BendingPhiVal[1] = getdPhi(recHitBendingPhi[SampleLayerCollection[1]][SampleLayerCollection[SampleLayerCollection.size()-1]], recHitBendingPhi[SampleLayerCollection[1]][SampleLayerCollection[1]]) * -1.;
            BendingPhiVal[2] = getdPhi(recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[SampleLayerCollection.size()-1]], recHitBendingPhi[SampleLayerCollection[0]][SampleLayerCollection[0]]) * -1.;
            BendingPhiVal[3] = recBendingPhiMax;
            // filter small reverse bending if needed
            if(fabs(recBendingPhiMax) < MaxBendingCut || fabs(BendingPhiVal[0]) < SegmentBendingCut0 || fabs(BendingPhiVal[1]) < SegmentBendingCut1 || fabs(BendingPhiVal[2]) < SegmentBendingCut2) {
                if(debug) cout << "block by BendingPhiTH." << endl;
                continue;
            }

            // for charge survey
            /*
            SimBendingPhiFilterSign = 0;
            for(int j = 0; j < 4; j++) {
                if(BendingPhiVal[j] != 0.)
                    SimBendingPhiValFilterSign[j] = (int)(BendingPhiVal[j]*(double)simTrackCharge/fabs(BendingPhiVal[j]));
                else
                    SimBendingPhiValFilterSign[j] = 0;
                SimBendingPhiFilterSign += SimBendingPhiValFilterSign[j];
            }
            //if(recBendingPhiMax*(double)simTrackCharge*BendingWiseCheck < 0.)
            if((double)SimBendingPhiFilterSign*BendingWiseCheck < 0. || (SimBendingPhiFilterSign == 0 && BendingWiseCheck != 0.))
                continue;
            */
            SimBendingPhiFilterSign = 0;
            for(int j = 0; j < 4; j++) {
                BendingPhiValSign[j] = 1;
                if((double)simTrackCharge * BendingPhiVal[j] <= 0.)
                    BendingPhiValSign[j] = 0;

                SimBendingPhiFilterSign += BendingPhiValSign[j];
            }
            if((double)(SimBendingPhiFilterSign-2.5)*BendingWiseCheck < 0. && BendingWiseCheck != 0.)
                continue;

            // check and correct reverse bending
            if(applyFilter == true) {
                BendingPhiFilterSign = 0;
                for(int j = 0; j < 4; j++) {
                    if(BendingPhiVal[j] != 0.)
                        BendingPhiValFilterSign[j] = (int)(BendingPhiVal[j]/fabs(BendingPhiVal[j]));
                    else
                        BendingPhiValFilterSign[j] = 0;
                    BendingPhiFilterSign += BendingPhiValFilterSign[j];
                }
                if(abs(BendingPhiFilterSign) < 2)
                    continue;

                //if((fabs(BendingPhiVal[0]) < 0.5*fabs(recBendingPhiMax) && BendingPhiValFilterSign[0]*BendingPhiFilterSign > 0.) || (fabs(BendingPhiVal[1]) > 0.9*fabs(recBendingPhiMax) && BendingPhiValFilterSign[1]*BendingPhiFilterSign > 0.) || (fabs(BendingPhiVal[2]) < 0.3*fabs(recBendingPhiMax) && BendingPhiValFilterSign[2]*BendingPhiFilterSign > 0.))
                    //continue;
                //if((BendingPhiValFilterSign[0]*BendingPhiFilterSign > 0. && BendingPhiValFilterSign[1]*BendingPhiFilterSign > 0. && BendingPhiValFilterSign[2]*BendingPhiFilterSign > 0.) && (fabs(BendingPhiVal[0]) < fabs(BendingPhiVal[1]) || fabs(BendingPhiVal[0]) < fabs(BendingPhiVal[2]) || fabs(BendingPhiVal[2]) < fabs(BendingPhiVal[1])))
                    //continue;

                //if(fabs(BendingPhiVal[0]/BendingPhiVal[3]) >= 0.99 && fabs(BendingPhiVal[3]) <= 0.01)
                    //continue;


                int badPattern = 0;
                if(BendingPhiValFilterSign[0]*BendingPhiFilterSign > 0. && BendingPhiValFilterSign[1]*BendingPhiFilterSign > 0. && fabs(BendingPhiVal[0]) < fabs(BendingPhiVal[1]))
                    badPattern++;
                if(BendingPhiValFilterSign[0]*BendingPhiFilterSign > 0. && BendingPhiValFilterSign[2]*BendingPhiFilterSign > 0. && fabs(BendingPhiVal[0]) < fabs(BendingPhiVal[2]))
                    badPattern++;
                if(BendingPhiValFilterSign[2]*BendingPhiFilterSign > 0. && BendingPhiValFilterSign[1]*BendingPhiFilterSign > 0. && fabs(BendingPhiVal[2]) < fabs(BendingPhiVal[1]))
                    badPattern++;
                //if(badPattern >= 3)
                    //continue;

                //for(int j = 0; j < 4; j++) {
                //    if((double)BendingPhiValFilterSign[j]*BendingPhiFilterSign <= 0.)
                //        BendingPhiVal[j] = 0.;
                //}
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

        if(ClusterSize[i] > 0 && ClusterSize[j] > 0 && ClusterSize[k] > 0 && ClusterSize[l] > 0 && ClusterSize[i] < 3 && ClusterSize[j] < 3 && ClusterSize[k] < 3 && ClusterSize[l] < 3) {
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
    recHitBendingPhiDSignHist->Fill(theRefPt, BendingPhiValSign[0]+BendingPhiValSign[1]+BendingPhiValSign[2]+BendingPhiValSign[3]);

    for(int i = 0; i < 3; i++) {
        if((applyFilter == true && (double)BendingPhiValFilterSign[i]*BendingPhiFilterSign > 0.) || applyFilter == false) {
            recHitBendingPhiSHist[i]->Fill(theRefPt, BendingPhiVal[i]);
            if(recBendingPhiMax != 0.)
                recHitBendingPhiDHist[i]->Fill(theRefPt, BendingPhiVal[i]/recBendingPhiMax);
        }
    }
    if((applyFilter == true && (double)BendingPhiValFilterSign[3]*BendingPhiFilterSign > 0.) || applyFilter == false)
        recHitBendingPhiSHist[3]->Fill(theRefPt, BendingPhiVal[3]);

    if(BendingPhiVal[0] != 0.) {
        recHitBendingPhiBHist[0]->Fill(theRefPt, BendingPhiVal[1]/BendingPhiVal[0]);
        recHitBendingPhiBHist[1]->Fill(theRefPt, BendingPhiVal[2]/BendingPhiVal[0]);
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
    /*
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
    */
    SaveName = FinalOutput + "recHitBendingDSignPhi" + Index + OutputFix;
    if(debug) cout << SaveName << endl;
    recHitBendingPhiDSignHist->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiDSignHist->GetXaxis()->CenterTitle(1);
    recHitBendingPhiDSignHist->GetYaxis()->SetTitle("DSign");
    recHitBendingPhiDSignHist->GetYaxis()->CenterTitle(1);
    recHitBendingPhiDSignHist->Draw("LEGO");
    BendingPhiCanvas->SaveAs(SaveName.c_str());
    
    for(int j = 0; j < 3; j++) {
        std::stringstream TempIndexJ;
        TempIndexJ << j;

        SaveName = FinalOutput + "recHitBendingPhiD" + TempIndexJ.str() + OutputFix;
        recHitBendingPhiDHist[j]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
        recHitBendingPhiDHist[j]->GetXaxis()->CenterTitle(1);
        recHitBendingPhiDHist[j]->GetYaxis()->SetTitle("Bending/Max");
        recHitBendingPhiDHist[j]->GetYaxis()->CenterTitle(1);
        recHitBendingPhiDHist[j]->Draw(theDrawOption.c_str());
        BendingPhiCanvas->SaveAs(SaveName.c_str());

        SaveName = FinalOutput + "recHitBendingPhiS" + TempIndexJ.str() + OutputFix;
        recHitBendingPhiSHist[j]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
        recHitBendingPhiSHist[j]->GetXaxis()->CenterTitle(1);
        recHitBendingPhiSHist[j]->GetYaxis()->SetTitle("Bending in Phi");
        recHitBendingPhiSHist[j]->GetYaxis()->CenterTitle(1);
        recHitBendingPhiSHist[j]->Draw(theDrawOption.c_str());
        BendingPhiCanvas->SaveAs(SaveName.c_str());
    }
    SaveName = FinalOutput + "recHitBendingPhiS3" + OutputFix;
    recHitBendingPhiSHist[3]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiSHist[3]->GetXaxis()->CenterTitle(1);
    recHitBendingPhiSHist[3]->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiSHist[3]->GetYaxis()->CenterTitle(1);
    recHitBendingPhiSHist[3]->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());

    SaveName = FinalOutput + "recHitBendingPhiB10" + OutputFix;
    recHitBendingPhiBHist[0]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiBHist[0]->GetXaxis()->CenterTitle(1);
    recHitBendingPhiBHist[0]->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiBHist[0]->GetYaxis()->CenterTitle(1);
    recHitBendingPhiBHist[0]->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());

    SaveName = FinalOutput + "recHitBendingPhiB20" + OutputFix;
    recHitBendingPhiBHist[1]->GetXaxis()->SetTitle("charge*simPt@Ref Gev");
    recHitBendingPhiBHist[1]->GetXaxis()->CenterTitle(1);
    recHitBendingPhiBHist[1]->GetYaxis()->SetTitle("Bending in Phi");
    recHitBendingPhiBHist[1]->GetYaxis()->CenterTitle(1);
    recHitBendingPhiBHist[1]->Draw(theDrawOption.c_str());
    BendingPhiCanvas->SaveAs(SaveName.c_str());
}

void RPCBendingAnalyzer::runSample(double ms0, double me0, double ss0, double se0, double ms1, double me1, double ss1, double se1, double ms2, double me2, double ss2, double se2, double ms3, double me3, double ss3, double se3) {

    int ms[4], me[4], ss[4],se[4];
    ms[0] = ms0; me[0] = me0; ss[0] = ss0; se[0] = se0;
    ms[1] = ms1; me[1] = me1; ss[1] = ss1; se[1] = se1;
    ms[2] = ms2; me[2] = me2; ss[2] = ss2; se[2] = se2;
    ms[3] = ms3; me[3] = me3; ss[3] = ss3; se[3] = se3;

    for(int Index = 0; Index < 4; Index++) {
        getEventRato(Index, 0);
        fitPtofPhi(Index, "landau", 0., ms[Index], me[Index], ss[Index], se[Index]);
    }
}

void RPCBendingAnalyzer::getEventRato(unsigned int IndexNumber, double theFitPtRange) {

    double FitPtRange = theFitPtRange;
    if(FitPtRange <= 0.)
        FitPtRange = WorkPtRange;

    TH2F *SampleHist = recHitBendingPhiSHist[IndexNumber];
    
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

    TH2F *SampleHist = recHitBendingPhiSHist[IndexNumber];

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
