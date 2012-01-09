#ifndef MyModule_BeamHaloAnalyzer_CSC2RPCBinder_H
#define MyModule_BeamHaloAnalyzer_CSC2RPCBinder_H

/** \class RPCSeedFinder
 *
 *   \author Haiyun.Teng - Peking University
 *
 *
 */

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>

using namespace std;
using namespace edm;

class CSC2RPCBinder {
    
    public:
        CSC2RPCBinder();
        ~CSC2RPCBinder();
        void initiate(edm::ESHandle<RPCGeometry>& RPCGeometryEH, edm::ESHandle<CSCGeometry>& CSCGeometryEH);
        void checkBinding(edm::ESHandle<RPCGeometry>& RPCGeometryEH, edm::ESHandle<CSCGeometry>& CSCGeometryEH);
        const std::vector<RPCDetId> centerRPCRolls(const CSCDetId& CSCId);
        const std::vector<RPCDetId> leftRPCRolls(const CSCDetId& CSCId);
        const std::vector<RPCDetId> rightRPCRolls(const CSCDetId& CSCId);
        const std::vector<RPCDetId> topRPCRolls(const CSCDetId& CSCId);
        const std::vector<RPCDetId> buttomRPCRolls(const CSCDetId& CSCId);
        const std::vector<RPCDetId> nearbyRPCRolls(const CSCDetId& CSCId);
        //const std::vector<RPCDetId> lefttopRPCRolls(const CSCDetId& CSCId);
        //const std::vector<RPCDetId> righttopRPCRolls(const CSCDetId& CSCId);
        //const std::vector<RPCDetId> leftbuttomRPCRolls(const CSCDetId& CSCId);
        //const std::vector<RPCDetId> rightbuttomRPCRolls(const CSCDetId& CSCId);
    private:
        bool initiated;
        // CSC strip number increases from left to right (1:n) when one looks from the small "inner" end toward the large "outer" end of the chamber
        std::map< CSCDetId, std::vector<RPCDetId> > centerRollswithCSCChamber;
        std::map< CSCDetId, std::vector<RPCDetId> > leftRollswithCSCChamber;
        std::map< CSCDetId, std::vector<RPCDetId> > rightRollswithCSCChamber;
        std::map< CSCDetId, std::vector<RPCDetId> > topRollswithCSCChamber;
        std::map< CSCDetId, std::vector<RPCDetId> > buttomRollswithCSCChamber;
        std::map< CSCDetId, std::vector<RPCDetId> > nearbyRollswithCSCChamber;
        //std::map< CSCDetId, std::vector<RPCDetId> > lefttopRollswithCSCChamber;
        //std::map< CSCDetId, std::vector<RPCDetId> > righttopRollswithCSCChamber;
        //std::map< CSCDetId, std::vector<RPCDetId> > leftbuttomRollswithCSCChamber;
        //std::map< CSCDetId, std::vector<RPCDetId> > rightbuttomRollswithCSCChamber;
};

#endif
