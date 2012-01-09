/*
 *  See header file for a description of this class.
 *
 *
 *  $Date: 2009/10/31 02:00:48 $
 *  \author Haiyun.Teng - Peking University
 *
 */

#include "MyModule/CSCSegment2RPC/src/CSC2RPCBinder.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <FWCore/Framework/interface/ESHandle.h>

#ifndef debug
#define debug 1
#endif


CSC2RPCBinder::CSC2RPCBinder() {
    initiated = false;
}

CSC2RPCBinder::~CSC2RPCBinder() {
}

void CSC2RPCBinder::initiate(edm::ESHandle<RPCGeometry>& RPCGeometryEH, edm::ESHandle<CSCGeometry>& CSCGeometryEH) {
    // Fill sub RPC map
    centerRollswithCSCChamber.clear();
    leftRollswithCSCChamber.clear();
    rightRollswithCSCChamber.clear();
    topRollswithCSCChamber.clear();
    buttomRollswithCSCChamber.clear();
    const std::vector<RPCRoll*> RPCRollPointerCollection= RPCGeometryEH->rolls();
    for(std::vector<RPCRoll*>::const_iterator RPCRollPointerIter = RPCRollPointerCollection.begin(); RPCRollPointerIter != RPCRollPointerCollection.end(); RPCRollPointerIter++) {
        RPCDetId RPCRollId = (*RPCRollPointerIter)->id();
        RPCGeomServ RPCServ(RPCRollId);
        int RPCRegion = RPCRollId.region();
        int RPCStation = RPCRollId.station();
        int RPCRing = RPCRollId.ring();
        int RPCChamber = RPCServ.segment();
        if(debug) cout << "RPC Roll RPCRegion " << RPCRegion << " RPCStation " << RPCStation << " RPCRing " << RPCRing << " RPCChamber " << RPCChamber << endl;
        int CSCRegion = RPCRegion;
        int CSCStation = RPCStation;
        int CSCRing = RPCRing;
        int CSCChamber = RPCChamber;
        int CSCLayer = 0; // Layer=0 for whole CSC chamber with 6 layers
        
        if(RPCRegion == 0)
            continue;
        // For RPC Region is +1/-1 for +Z/-Z, while for CSC endcap is 1/2 for +Z/-Z
        if(RPCRegion == -1)
            CSCRegion = 2;
        if(RPCStation == 4)
            continue;
        if((RPCStation == 2 || RPCStation == 3) && RPCRing == 3)
            CSCRing = 2;
        
        // For center CSC
        CSCDetId CSCIndex(CSCRegion, CSCStation, CSCRing, CSCChamber, CSCLayer);
        std::vector<RPCDetId> tempRPCRolls;
        if(centerRollswithCSCChamber.find(CSCIndex) != centerRollswithCSCChamber.end())
            tempRPCRolls = centerRollswithCSCChamber[CSCIndex];
        tempRPCRolls.push_back(RPCRollId);
        centerRollswithCSCChamber[CSCIndex] = tempRPCRolls;

        // For aside CSC
        if(RPCStation == 1 && RPCRing == 2) {
            int leftChamber = (CSCChamber+RPCRegion+35)%36+1;
            int rightChamber = (CSCChamber-RPCRegion+35)%36+1;
            CSCDetId CSCIndex0(CSCRegion, CSCStation, CSCRing-1, CSCChamber, CSCLayer);
            CSCDetId CSCIndex1(CSCRegion, CSCStation, CSCRing+1, CSCChamber, CSCLayer);
            CSCDetId CSCIndex2(CSCRegion, CSCStation, CSCRing, leftChamber, CSCLayer);
            CSCDetId CSCIndex3(CSCRegion, CSCStation, CSCRing, rightChamber, CSCLayer);
            
            if(debug) cout << "leftChamber is " << leftChamber << ", rightChamber is " << rightChamber << endl;

            std::vector<RPCDetId> RPCRolls_temp0;
            if(topRollswithCSCChamber.find(CSCIndex0) != topRollswithCSCChamber.end())
                RPCRolls_temp0 = topRollswithCSCChamber[CSCIndex0];
            RPCRolls_temp0.push_back(RPCRollId);
            topRollswithCSCChamber[CSCIndex0] = RPCRolls_temp0;

            std::vector<RPCDetId> RPCRolls_temp1;
            if(buttomRollswithCSCChamber.find(CSCIndex1) != buttomRollswithCSCChamber.end())
                RPCRolls_temp1 = buttomRollswithCSCChamber[CSCIndex1];
            RPCRolls_temp1.push_back(RPCRollId);
            buttomRollswithCSCChamber[CSCIndex1] = RPCRolls_temp1;

            std::vector<RPCDetId> RPCRolls_temp2;
            if(leftRollswithCSCChamber.find(CSCIndex2) != leftRollswithCSCChamber.end())
                RPCRolls_temp2 = leftRollswithCSCChamber[CSCIndex2];
            RPCRolls_temp2.push_back(RPCRollId);
            leftRollswithCSCChamber[CSCIndex2] = RPCRolls_temp2;

            std::vector<RPCDetId> RPCRolls_temp3;
            if(rightRollswithCSCChamber.find(CSCIndex3) != rightRollswithCSCChamber.end())
                RPCRolls_temp3 = rightRollswithCSCChamber[CSCIndex3];
            RPCRolls_temp3.push_back(RPCRollId);
            rightRollswithCSCChamber[CSCIndex3] = RPCRolls_temp3;
        }
        
        if(RPCStation == 1 && RPCRing == 3) {
            int leftChamber = (CSCChamber+RPCRegion+35)%36+1;
            int rightChamber = (CSCChamber-RPCRegion+35)%36+1;
            CSCDetId CSCIndex0(CSCRegion, CSCStation, CSCRing-1, CSCChamber, CSCLayer);
            CSCDetId CSCIndex2(CSCRegion, CSCStation, CSCRing, leftChamber, CSCLayer);
            CSCDetId CSCIndex3(CSCRegion, CSCStation, CSCRing, rightChamber, CSCLayer);

            if(debug) cout << "leftChamber is " << leftChamber << ", rightChamber is " << rightChamber << endl;

            std::vector<RPCDetId> RPCRolls_temp0;
            if(topRollswithCSCChamber.find(CSCIndex0) != topRollswithCSCChamber.end())
                RPCRolls_temp0 = topRollswithCSCChamber[CSCIndex0];
            RPCRolls_temp0.push_back(RPCRollId);
            topRollswithCSCChamber[CSCIndex0] = RPCRolls_temp0;

            std::vector<RPCDetId> RPCRolls_temp2;
            if(leftRollswithCSCChamber.find(CSCIndex2) != leftRollswithCSCChamber.end())
                RPCRolls_temp2 = leftRollswithCSCChamber[CSCIndex2];
            RPCRolls_temp2.push_back(RPCRollId);
            leftRollswithCSCChamber[CSCIndex2] = RPCRolls_temp2;

            std::vector<RPCDetId> RPCRolls_temp3;
            if(rightRollswithCSCChamber.find(CSCIndex3) != rightRollswithCSCChamber.end())
                RPCRolls_temp3 = rightRollswithCSCChamber[CSCIndex3];
            RPCRolls_temp3.push_back(RPCRollId);
            rightRollswithCSCChamber[CSCIndex3] = RPCRolls_temp3;
        }

        if(RPCStation != 1 && RPCRing == 2) {
            int leftChamber = (CSCChamber+RPCRegion*((CSCStation<=2)?1:-1)+35)%36+1;
            int rightChamber = (CSCChamber-RPCRegion*((CSCStation<=2)?1:-1)+35)%36+1;
            CSCDetId CSCIndex0(CSCRegion, CSCStation, CSCRing-1, (CSCChamber+1)/2, CSCLayer);
            CSCDetId CSCIndex2(CSCRegion, CSCStation, CSCRing, leftChamber, CSCLayer);
            CSCDetId CSCIndex3(CSCRegion, CSCStation, CSCRing, rightChamber, CSCLayer);

            if(debug) cout << "leftChamber is " << leftChamber << ", rightChamber is " << rightChamber << endl;
            
            std::vector<RPCDetId> RPCRolls_temp0;
            if(topRollswithCSCChamber.find(CSCIndex0) != topRollswithCSCChamber.end())
                RPCRolls_temp0 = topRollswithCSCChamber[CSCIndex0];
            RPCRolls_temp0.push_back(RPCRollId);
            topRollswithCSCChamber[CSCIndex0] = RPCRolls_temp0;

            std::vector<RPCDetId> RPCRolls_temp2;
            if(leftRollswithCSCChamber.find(CSCIndex2) != leftRollswithCSCChamber.end())
                RPCRolls_temp2 = leftRollswithCSCChamber[CSCIndex2];
            RPCRolls_temp2.push_back(RPCRollId);
            leftRollswithCSCChamber[CSCIndex2] = RPCRolls_temp2;

            std::vector<RPCDetId> RPCRolls_temp3;
            if(rightRollswithCSCChamber.find(CSCIndex3) != rightRollswithCSCChamber.end())
                RPCRolls_temp3 = rightRollswithCSCChamber[CSCIndex3];
            RPCRolls_temp3.push_back(RPCRollId);
            rightRollswithCSCChamber[CSCIndex3] = RPCRolls_temp3;
        }
        
        if(RPCStation != 1 && RPCRing == 3) {
            int leftChamber = (CSCChamber+RPCRegion*((CSCStation<=2)?1:-1)+35)%36+1;
            int rightChamber = (CSCChamber-RPCRegion*((CSCStation<=2)?1:-1)+35)%36+1;
            CSCDetId CSCIndex2(CSCRegion, CSCStation, CSCRing, leftChamber, CSCLayer);
            CSCDetId CSCIndex3(CSCRegion, CSCStation, CSCRing, rightChamber, CSCLayer);

            if(debug) cout << "leftChamber is " << leftChamber << ", rightChamber is " << rightChamber << endl;

            std::vector<RPCDetId> RPCRolls_temp2;
            if(leftRollswithCSCChamber.find(CSCIndex2) != leftRollswithCSCChamber.end())
                RPCRolls_temp2 = leftRollswithCSCChamber[CSCIndex2];
            RPCRolls_temp2.push_back(RPCRollId);
            leftRollswithCSCChamber[CSCIndex2] = RPCRolls_temp2;

            std::vector<RPCDetId> RPCRolls_temp3;
            if(rightRollswithCSCChamber.find(CSCIndex3) != rightRollswithCSCChamber.end())
                RPCRolls_temp3 = rightRollswithCSCChamber[CSCIndex3];
            RPCRolls_temp3.push_back(RPCRollId);
            rightRollswithCSCChamber[CSCIndex3] = RPCRolls_temp3;
        }
    }
    
    // Fill all nearby RPC map
    nearbyRollswithCSCChamber.clear();
    const std::vector<DetId> CSCDetIds = CSCGeometryEH->detIds();
    for(std::vector<DetId>::const_iterator CSCDetIdIter = CSCDetIds.begin(); CSCDetIdIter != CSCDetIds.end(); CSCDetIdIter++) {
        CSCDetId CSCIndex = CSCDetIdIter->rawId();
        // Center Rolls are always the 1st to be check
        nearbyRollswithCSCChamber[CSCIndex] = centerRollswithCSCChamber[CSCIndex];
        for(std::vector<RPCDetId>::const_iterator RPCDetIdIter = leftRollswithCSCChamber[CSCIndex].begin(); RPCDetIdIter != leftRollswithCSCChamber[CSCIndex].end(); RPCDetIdIter++)
            nearbyRollswithCSCChamber[CSCIndex].push_back(*RPCDetIdIter);
        for(std::vector<RPCDetId>::const_iterator RPCDetIdIter = rightRollswithCSCChamber[CSCIndex].begin(); RPCDetIdIter != rightRollswithCSCChamber[CSCIndex].end(); RPCDetIdIter++)
            nearbyRollswithCSCChamber[CSCIndex].push_back(*RPCDetIdIter);
        // we have check Y axis need to be constrait to center part, so we don't need the top/buttom binding
        /*
        for(std::vector<RPCDetId>::const_iterator RPCDetIdIter = topRollswithCSCChamber[CSCIndex].begin(); RPCDetIdIter != topRollswithCSCChamber[CSCIndex].end(); RPCDetIdIter++)
            nearbyRollswithCSCChamber[CSCIndex].push_back(*RPCDetIdIter);
        for(std::vector<RPCDetId>::const_iterator RPCDetIdIter = buttomRollswithCSCChamber[CSCIndex].begin(); RPCDetIdIter != buttomRollswithCSCChamber[CSCIndex].end(); RPCDetIdIter++)
            nearbyRollswithCSCChamber[CSCIndex].push_back(*RPCDetIdIter);
        */
    }
}

const std::vector<RPCDetId> CSC2RPCBinder::centerRPCRolls(const CSCDetId& CSCId) {
    return centerRollswithCSCChamber[CSCId];
}
        
const std::vector<RPCDetId> CSC2RPCBinder::leftRPCRolls(const CSCDetId& CSCId) {
    return leftRollswithCSCChamber[CSCId];
}
        
const std::vector<RPCDetId> CSC2RPCBinder::rightRPCRolls(const CSCDetId& CSCId) {
    return rightRollswithCSCChamber[CSCId];
}
        
const std::vector<RPCDetId> CSC2RPCBinder::topRPCRolls(const CSCDetId& CSCId) {
    return topRollswithCSCChamber[CSCId];
}
        
const std::vector<RPCDetId> CSC2RPCBinder::buttomRPCRolls(const CSCDetId& CSCId) {
    return buttomRollswithCSCChamber[CSCId];
}

const std::vector<RPCDetId> CSC2RPCBinder::nearbyRPCRolls(const CSCDetId& CSCId) {
    return nearbyRollswithCSCChamber[CSCId];
}

void CSC2RPCBinder::checkBinding(edm::ESHandle<RPCGeometry>& RPCGeometryEH, edm::ESHandle<CSCGeometry>& CSCGeometryEH) {    
    // Check the binding for center RPC chamber
    const std::vector<CSCChamber*> cscChambers = CSCGeometryEH->chambers();
    for(std::vector<CSCChamber*>::const_iterator CSCChamberIter = cscChambers.begin(); CSCChamberIter != cscChambers.end(); CSCChamberIter++) {
        CSCDetId cscChamber = (*CSCChamberIter)->id();
        int Endcap = cscChamber.endcap();
        int Station = cscChamber.station();
        int Ring = cscChamber.ring();
        int Chamber = cscChamber.chamber();
        int Layer = cscChamber.layer();
        int PhivsStripDirection = ((Station <= 2) ? 1 : -1) * ((Endcap == 1) ? 1 : -1);
        if(debug) cout << "Checking endcap CSC chamber: " << cscChamber.rawId() << ", Endcap: " << Endcap << ", Station: " << Station << ", Ring: " << Ring << ", Chamber: " << Chamber << ", layer: " << Layer << endl;
        //const Bounds& cscBounds = (*CSCChamberIter)->surface().bounds();
        const CSCLayer* cscLayer = (*CSCChamberIter)->layer(1);
        const CSCStripTopology* cscTop = cscLayer->geometry()->topology();
        int nCSCStrip0 = cscLayer->geometry()->numberOfStrips();
        int nCSCStrip1 = cscTop->nstrips();
        if(debug) cout << "CSC layer strips from layer geometry: " << nCSCStrip0 << ", from layer topology: " << nCSCStrip1 << endl;
        GlobalPoint CSCEdge0 = cscLayer->toGlobal(cscTop->localPosition(0));
        GlobalPoint CSCEdge1 = cscLayer->toGlobal(cscTop->localPosition(nCSCStrip1));
        if(debug) cout << "CSC Phi range from " << CSCEdge0.phi() << " to " << CSCEdge1.phi() << endl;
        if(debug) cout << "CSC Phi degree range from " << CSCEdge0.phi().degrees() << " to " << CSCEdge1.phi().degrees() << endl;
        if(centerRollswithCSCChamber.find(cscChamber) != centerRollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = centerRollswithCSCChamber[cscChamber];
            if(debug) cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                const RPCRoll* rpcRoll = RPCGeometryEH->roll(*RPCRollIter);
                int nRPCStrip = rpcRoll->nstrips();
                double RPCPitch = rpcRoll->pitch();
                GlobalPoint RPCEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0));
                GlobalPoint RPCEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip));
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << RPCEdge0.phi() << " to " << RPCEdge1.phi() << endl;
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << RPCEdge0.phi().degrees() << " to " << RPCEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge0).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge1).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge1).y(), 0);
                if(debug) cout << "RPC roll nstrips: " << nRPCStrip << ", pitch: " << RPCPitch << ", roll size: " << RPCPitch*nRPCStrip << endl;
                //if(cscBounds.inside(inside0) && cscBounds.inside(inside1))
                if((CSCEdge0.phi()-RPCEdge0.phi()).value()*(CSCEdge1.phi()-RPCEdge1.phi()).value() < 0) {
                    if(debug) cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
                else {
                    if(debug) cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
            }
        }
        else {
            if(debug) cout << "Could not find the center binding RPC roll" << endl;
        }

        // Check the binding for left RPC chamber
        if(centerRollswithCSCChamber.find(cscChamber) != centerRollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = leftRollswithCSCChamber[cscChamber];
            if(debug) cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                const RPCRoll* rpcRoll = RPCGeometryEH->roll(*RPCRollIter);
                int nRPCStrip = rpcRoll->nstrips();
                GlobalPoint RPCEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0));
                GlobalPoint RPCEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip));
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << RPCEdge0.phi() << " to " << RPCEdge1.phi() << endl;
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << RPCEdge0.phi().degrees() << " to " << RPCEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge0).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge1).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge1).y(), 0);
                if(((CSCEdge0.phi()-RPCEdge0.phi()).value()*PhivsStripDirection) > 0 && ((CSCEdge1.phi()-RPCEdge1.phi()).value()*PhivsStripDirection) > 0) {
                    if(debug) cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
                else {
                    if(debug) cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
            }
        }
        else {
            if(debug) cout << "Could not find the binding left RPC roll" << endl;
        }

        // Check the binding for right RPC chamber
        if(centerRollswithCSCChamber.find(cscChamber) != centerRollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = rightRollswithCSCChamber[cscChamber];
            if(debug) cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                const RPCRoll* rpcRoll = RPCGeometryEH->roll(*RPCRollIter);
                int nRPCStrip = rpcRoll->nstrips();
                GlobalPoint RPCEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0));
                GlobalPoint RPCEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip));
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << RPCEdge0.phi() << " to " << RPCEdge1.phi() << endl;
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << RPCEdge0.phi().degrees() << " to " << RPCEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge0).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge1).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge1).y(), 0);
                if(((CSCEdge0.phi()-RPCEdge0.phi()).value()*PhivsStripDirection) < 0 && ((CSCEdge1.phi()-RPCEdge1.phi()).value()*PhivsStripDirection) < 0) {
                    if(debug) cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
                else {
                    if(debug) cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
            } 
        } 
        else {
            if(debug) cout << "Could not find the binding right RPC roll" << endl;
        }

        // Check the binding for top RPC chamber
        if(centerRollswithCSCChamber.find(cscChamber) != centerRollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = topRollswithCSCChamber[cscChamber];
            if(debug) cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                int theCSCRing = rRing;
                if(rStation != 1 && rRing == 3)
                    theCSCRing = 2;
                const RPCRoll* rpcRoll = RPCGeometryEH->roll(*RPCRollIter);
                int nRPCStrip = rpcRoll->nstrips();
                GlobalPoint RPCEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0));
                GlobalPoint RPCEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip));
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << RPCEdge0.phi() << " to " << RPCEdge1.phi() << endl;
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << RPCEdge0.phi().degrees() << " to " << RPCEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge0).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge1).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge1).y(), 0);
                if((CSCEdge0.phi()-RPCEdge0.phi()).value()*(CSCEdge1.phi()-RPCEdge1.phi()).value() < 0 || (theCSCRing - Ring) == 1) {
                    if(debug) cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
                else {
                    if(debug) cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
            } 
        } 
        else {
            if(debug) cout << "Could not find the binding top RPC roll" << endl;
        }

        // Check the binding for buttom RPC chamber
        if(centerRollswithCSCChamber.find(cscChamber) != centerRollswithCSCChamber.end()) {
            std::vector<RPCDetId> RPCRolls;
            RPCRolls = buttomRollswithCSCChamber[cscChamber];
            if(debug) cout <<"Bind to " << RPCRolls.size() << " RPC Rolls" << endl;
            for(std::vector<RPCDetId>::const_iterator RPCRollIter = RPCRolls.begin(); RPCRollIter != RPCRolls.end(); RPCRollIter++) {
                RPCGeomServ rpcSrv(*RPCRollIter);
                int rEndcap = RPCRollIter->region();
                int rStation = RPCRollIter->station();
                int rRing = RPCRollIter->ring();
                int rSegment = rpcSrv.segment();
                int rRoll = RPCRollIter->roll();
                int theCSCRing = rRing;
                if(rStation != 1 && rRing == 3)
                    theCSCRing = 2;
                const RPCRoll* rpcRoll = RPCGeometryEH->roll(*RPCRollIter); 
                int nRPCStrip = rpcRoll->nstrips(); 
                GlobalPoint RPCEdge0 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(0)); 
                GlobalPoint RPCEdge1 = rpcRoll->toGlobal(rpcRoll->centreOfStrip(nRPCStrip)); 
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi range from " << RPCEdge0.phi() << " to " << RPCEdge1.phi() << endl;
                if(debug) cout << "RPC roll: " << RPCRollIter->rawId() <<" Phi degree range from " << RPCEdge0.phi().degrees() << " to " << RPCEdge1.phi().degrees() << endl;
                LocalPoint inside0 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge0).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge0).y(), 0);
                LocalPoint inside1 = LocalPoint((*CSCChamberIter)->surface().toLocal(RPCEdge1).x(), (*CSCChamberIter)->surface().toLocal(RPCEdge1).y(), 0);
                if((CSCEdge0.phi()-RPCEdge0.phi()).value()*(CSCEdge1.phi()-RPCEdge1.phi()).value() < 0 || (theCSCRing - Ring) == -1) {
                    if(debug) cout << "Well binding for RPC Roll: " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
                else {
                    if(debug) cout << "Wrong binding for RPC Roll " << RPCRollIter->rawId() << ", Endcap: " << rEndcap << ", Station: " << rStation << ", Ring: " << rRing << ", Segment: " << rSegment << ", Roll: " << rRoll << endl;
                }
            } 
        } 
        else {
            if(debug) cout << "Could not find the binding buttom RPC roll" << endl;
        }
    }
}
