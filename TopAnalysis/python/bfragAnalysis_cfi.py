import FWCore.ParameterSet.Config as cms

bfragAnalysis = cms.EDAnalyzer("FragmentationAnalyzer",
                               hadronList   = cms.vint32( 511 ,521,531 ,5122,413 ,423 ,431,411,10411,10421,10413,10423,20413,20423,415,425,10431,433,10433,20433,435),
                               DssList      = cms.vint32( 411,10411,10421,10413,10423,20413,20423,415,425,10431,433,10433,20433,435),
                               hadronBSF    = cms.vdouble(0.941, 0.942, 1.088, 1.87),
                               hadronBunc   = cms.vdouble(0.014, 0.014, 0.052,0.26),
                               hadronDzSF   = cms.vdouble(0.780, 1.020, 0.857, 2.90),
                               hadronDzuuc  = cms.vdouble(0.026, 0.027, 0.054, 0.24),
                               hadronDzb    = cms.vdouble(47.4,74 ,1.4 ,1.1 ,0.5 ,74.7,0.4),
                               hadronUncDzb = cms.vdouble(2.8 ,4  ,0.13,1.1 ,67.7,0.9 ,0.4),
                               hadronDz     = cms.vdouble(8.1 ,8.6,1.9 ,10.9,67.7,64.7,0.9,0.4),
                               hadronUncDz  = cms.vdouble(1.5 ,0.7,0.5 ,1.1 ,0.5 ,0.9 ,0.9,0.4),
                               hadronUncJPsi= cms.vdouble(0.0028, 0.0034, 0.0021, 0.0003, 0, 0, 0, 0),
                               numEntries = cms.vint32(500000),
                               GENInfo      = cms.InputTag("generator")
                               )
