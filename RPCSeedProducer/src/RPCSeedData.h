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


#ifndef BarrelDoubleSegmentCode1
#define BarrelDoubleSegmentCode1 47 // 101111
#endif

#ifndef BarrelDoubleSegmentCode2
#define BarrelDoubleSegmentCode2 31 // 011111
#endif

#ifndef BarrelDoubleSegmentOptionalCode
#define BarrelDoubleSegmentOptionalCode 48 // 110000
#endif

#ifndef BarrelSingleSegmentCode1
#define BarrelSingleSegmentCode1 51 // 110011
#endif

#ifndef BarrelSingleSegmentCode2
#define BarrelSingleSegmentCode2 60 // 111100
#endif

#ifndef BarrelSingleSegmentOptionalCode1
#define BarrelSingleSegmentOptionalCode1 3 // 000011
#endif

#ifndef BarrelSingleSegmentOptionalCode2
#define BarrelSingleSegmentOptionalCode2 12 // 001100
#endif

#ifndef BarrelSingleSegmentOptionalCode3
#define BarrelSingleSegmentOptionalCode3 16 // 010000
#endif

#ifndef BarrelSingleSegmentOptionalCode4
#define BarrelSingleSegmentOptionalCode4 32 // 100000
#endif

#ifndef EndcapSingleSegmentCode1
#define EndcapSingleSegmentCode1 15 // 1111
#endif

#ifndef EndcapSingleSegmentOptionalCode
#define EndcapSingleSegmentOptionalCode 24 // 11000
#endif

#ifndef UpperLimitMeanPt2
#define UpperLimitMeanPt2 50.0
#endif

#ifndef UpperLimitSigmaPt2
#define UpperLimitSigmaPt2 20.0
#endif

#ifndef UpperLimitMeanPt1
#define UpperLimitMeanPt1 20.0
#endif

#ifndef UpperLimitSigmaPt1 
#define UpperLimitSigmaPt1 5.0
#endif

#ifndef LowerLimitMeanPt
#define LowerLimitMeanPt 4.0
#endif

#ifndef LowerLimitSigmaPt
#define LowerLimitSigmaPt 1.0
#endif

#endif
