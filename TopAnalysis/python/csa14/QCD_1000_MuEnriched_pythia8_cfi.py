import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/109C58A8-B705-E411-8716-00261894386F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/443C36E6-B805-E411-9AE2-0025905A60AA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/449220CB-0B06-E411-BB24-0025905A612C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/52226338-B405-E411-866B-00261894388F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5650BBDE-B005-E411-825D-0026189438DD.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/66460EAB-AE05-E411-88C4-0025905A6066.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/88C61653-C005-E411-A437-0025905A6082.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/9A05CCC4-0B06-E411-BADB-002590596498.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/9A39A592-B305-E411-B6FE-002618943867.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B6F76AE6-B005-E411-A609-0026189438F8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C4CBCE49-0C06-E411-BE47-003048FFD760.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E6235CE9-B005-E411-B26E-002618943982.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E678EFF7-C305-E411-A7B5-002590596498.root' ] );


secFiles.extend( [
               ] )
