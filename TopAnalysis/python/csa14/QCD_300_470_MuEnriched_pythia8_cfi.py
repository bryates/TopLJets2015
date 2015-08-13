import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/02138B8A-FD05-E411-ACD3-001F2965D276.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0274BA9C-FB05-E411-A1FC-001F2965648A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/081CDCEA-FB05-E411-A76D-002590D601B8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/14C71BFE-FB05-E411-888F-D48564599C64.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1CEC8E34-FD05-E411-8043-D48564593FA8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/26B17062-FB05-E411-A5F4-00259073E520.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4676D284-FB05-E411-B126-20CF3027A5F4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4892D6D5-FB05-E411-B6D2-D48564593FA8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/48F514C1-FB05-E411-815B-00259073E39C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4CCD502B-FC05-E411-B109-002590D8C746.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/56C518DF-FB05-E411-AE43-D48564594FB4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/76ACE984-FB05-E411-A6C8-002590207C28.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/785FDE8A-FB05-E411-B2FD-20CF307C98DC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/86832FBF-FD05-E411-A568-002590207C28.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/882EF341-FB05-E411-BB33-001F2965348A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AC522559-FB05-E411-8DA0-00259073E3B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B0C1CBD5-FB05-E411-A829-001F29659BAE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B8CD28BF-FB05-E411-9421-0025907277CE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BAEB7BF6-FB05-E411-A95E-001F2965D276.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/CA0E64BD-FB05-E411-B3D4-001F2965444E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D06B8E75-FC05-E411-B108-00259073E544.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D69AA296-FB05-E411-A097-D48564599CEE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E41EAE70-FC05-E411-A177-D48564594F36.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E8CB8EAB-FD05-E411-9DF0-20CF307C98DC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/FA3D383A-FB05-E411-9A2C-F4CE46B27A1A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-300to470_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/FE55D737-FC05-E411-8BDF-00259073E468.root' ] );


secFiles.extend( [
               ] )

