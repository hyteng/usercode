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

// RPCLayer: 5/4/3/2/1/0
#ifndef BarrelPatternCode1
#define BarrelPatternCode1 15 // 001111=RPCLayer0123
#endif

#ifndef BarrelPatternCode2
#define BarrelPatternCode2 23 // 010111=RPCLayer0124
#endif

#ifndef BarrelPatternCode3
#define BarrelPatternCode3 39 // 100111=RPCLayer0125
#endif

#ifndef BarrelPatternCode4
#define BarrelPatternCode4 27 // 011011=RPCLayer0134
#endif

#ifndef BarrelPatternCode5
#define BarrelPatternCode5 43 // 101011=RPCLayer0135
#endif

#ifndef BarrelPatternCode6
#define BarrelPatternCode6 51 // 110011=RPCLayer0145
#endif

#ifndef BarrelPatternCode7
#define BarrelPatternCode7 29 // 011101=RPCLayer0234
#endif

#ifndef BarrelPatternCode8
#define BarrelPatternCode8 45 // 101101=RPCLayer0235
#endif

#ifndef BarrelPatternCode9
#define BarrelPatternCode9 53 // 110101=RPCLayer0245
#endif

#ifndef BarrelPatternCode10
#define BarrelPatternCode10 57 // 111001=RPCLayer0345
#endif

#ifndef BarrelPatternCode11
#define BarrelPatternCode11 30 // 011110=RPCLayer1234
#endif

#ifndef BarrelPatternCode12
#define BarrelPatternCode12 46 // 101110=RPCLayer1235
#endif

#ifndef BarrelPatternCode13
#define BarrelPatternCode13 54 // 110110=RPCLayer1245
#endif

#ifndef BarrelPatternCode14
#define BarrelPatternCode14 58 // 111010=RPCLayer1345
#endif

#ifndef BarrelPatternCode15
#define BarrelPatternCode15 60 // 111100=RPCLayer2345
#endif


#ifndef UpperLimitMeanPt2
#define UpperLimitMeanPt2 50.0
#endif

#ifndef UpperLimitSigmaPt2
#define UpperLimitSigmaPt2 10.0
#endif

#ifndef UpperLimitMeanPt1
#define UpperLimitMeanPt1 30.0
#endif

#ifndef UpperLimitSigmaPt1 
#define UpperLimitSigmaPt1 6.0
#endif

#ifndef LowerLimitMeanPt
#define LowerLimitMeanPt 4.0
#endif

#ifndef LowerLimitSigmaPt
#define LowerLimitSigmaPt 1.0
#endif

#ifndef PI
#define PI 3.1415926535
#endif

#endif
