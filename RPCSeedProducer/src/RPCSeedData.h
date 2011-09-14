#ifndef MyModule_RPCSeedProducer_RPCData_H
#define MyModule_RPCSeedProducer_RPCData_H

// in future RPC endcap may have 5 layers, while implant the double layer in RE2
#ifndef BarrelLayerNumber
#define BarrelLayerNumber 6
#endif

#ifndef EachEndcapLayerNumber
#define EachEndcapLayerNumber 3
#endif

#ifndef RPCLayerNumber
#define RPCLayerNumber 12
#endif

#ifndef NegativeEndcapLayerShift
#define NegativeEndcapLayerShift 5
#endif

#ifndef PositiveEndcapLayerShift
#define PositiveEndcapLayerShift 8
#endif


#ifndef BarrelDoubleSegmentCode
#define BarrelDoubleSegmentCode 15 // 001111
#endif
 
#ifndef BarrelDoubleSegmentOptionalCode
#define BarrelDoubleSegmentOptionalCode 48 // 110000
#endif

#ifndef BarrelSingleSegmentCode
#define BarrelSingleSegmentCode 11 // 001011
#endif

#ifndef BarrelSingleSegmentOptionalCode
#define BarrelSingleSegmentOptionalCode 48 // 110000
#endif

#ifndef EndcapSingleSegmentCode
#define EndcapSingleSegmentCode 6 // 00110
#endif

#ifndef EndcapSingleSegmentOptionalCode
#define EndcapSingleSegmentOptionalCode 24 // 11000
#endif

#endif
