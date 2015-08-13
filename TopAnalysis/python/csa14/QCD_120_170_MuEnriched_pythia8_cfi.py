import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/02904D1E-D622-E411-AB19-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/066B910F-D622-E411-BD66-0025905B85AE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/1432F120-D622-E411-B459-0025905A6080.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/28447D41-D622-E411-A74B-0025905A48EC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/2A7A7741-D622-E411-900D-0025905A48EC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/32CD4E1E-D622-E411-980E-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/387B450A-D622-E411-8894-0025905A6066.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/3CAB481E-D622-E411-A3CE-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/4080E70A-D622-E411-AA84-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/4E1E7841-D622-E411-8BA2-0025905A48EC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/5E1CEA20-D622-E411-872E-0025905A6080.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/5E3C220F-D622-E411-9D98-0025905B85AE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/62A54308-D622-E411-B67B-0025905B857C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/76AFF00A-D622-E411-AE66-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/7CC17941-D622-E411-A68C-0025905A48EC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/7EF1F10A-D622-E411-B7EF-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/84D9EA0A-D622-E411-BC2F-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/8AA9381E-D622-E411-AFEF-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/8CF25D08-D622-E411-9AFD-0025905B85B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/94142A21-D622-E411-937D-0025905A6080.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/B0B6C111-D622-E411-BDA9-0025905A60D6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/B6871906-D622-E411-8A55-0025905B855E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/D247E80A-D622-E411-A532-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/D6057C41-D622-E411-9813-0025905A48EC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/F00DEE20-D622-E411-9C2F-0025905A6080.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/F2D2E70A-D622-E411-8206-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/FAA88241-D622-E411-B754-0025905A48EC.root' ] );


secFiles.extend( [
               ] )


