import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/24357426-2709-E411-BDA0-18A905709A8A.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4CBB8736-2809-E411-A518-18A905704650.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5C7FD735-2809-E411-9219-18A9057049CC.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/70420636-2709-E411-9D4C-FA163EE95DCF.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/747E6424-2709-E411-88FB-18A905706BA8.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AAD39935-2809-E411-9F50-18A905708AA2.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/ECFD2C36-2809-E411-A633-18A90570ABE0.root' ] );


secFiles.extend( [
               ] )
