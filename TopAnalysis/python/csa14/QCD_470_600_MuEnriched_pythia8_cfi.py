import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/2ABDAFD0-1119-E411-8B0C-FA163E95E612.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/2E9012AB-1119-E411-BCDE-FA163E4AE49D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/5A4A923E-2219-E411-85B5-02163E00E8EB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/666588D2-1119-E411-A261-FA163EEA7264.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/82AA4F8F-2219-E411-9479-02163E00A123.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/94B00B89-2219-E411-99B2-02163E00EA7C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/9C51373C-2219-E411-98AC-02163E008DA1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/AA93F325-2219-E411-98B4-02163E00EABD.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/BE415DC4-1119-E411-A6A7-FA163E08C541.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/D868C089-1119-E411-B5C0-FA163E122695.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/D88EA30B-2219-E411-817C-02163E008DA1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/E6421622-2219-E411-8708-02163E00F1E5.root',
       '/store/mc/Spring14miniaod/QCD_Pt-470to600_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/30000/FE720B1C-2219-E411-B72E-02163E00ECD0.root' ] );


secFiles.extend( [
               ] )

