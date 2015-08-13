import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/00512A93-5619-E411-8C8A-001E6739BC29.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/1E4DDBD3-5919-E411-B2AF-001E673C84B9.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/2649B04F-5A19-E411-AB9A-001E6739CEB1.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/2CA8E3D9-5F19-E411-8095-9CB65404FBA0.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/54C3FB6E-6019-E411-BAAB-9CB65404ED04.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/5A6F4E35-5919-E411-92D6-001E67496A6C.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/6491397A-4D19-E411-ACC7-001E673D23F9.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/68FE1478-5819-E411-9BF1-001E67496A6C.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/6E934C11-5619-E411-9246-001E673C7EF4.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/70DCF56A-5819-E411-96B2-001E6739C8B9.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/74DEF2D8-5019-E411-95B3-38EAA7A6DAC8.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/8817DB6A-4D19-E411-B948-001E6739A781.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/983B1DA2-5919-E411-8DB6-001E6739B019.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/A422570A-5619-E411-B3D8-001E6739BC29.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/ACF8AFB1-5919-E411-A8BA-38EAA7A6D65C.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/BC1EC755-5919-E411-B7AD-38EAA7A6DAA0.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/BEDE4B01-A819-E411-8AFE-001E673C8767.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/CEF3A826-5119-E411-A892-001E673C84B9.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/D4F6BD1D-5A19-E411-9543-9CB65404EEF0.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/D6E19880-4D19-E411-89E9-001E6739B019.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/DC3F9094-5619-E411-B0EE-38EAA7A6D65C.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/E405AEFA-4F19-E411-BBFC-38EAA7A6DCF0.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/E4EDF628-5019-E411-9BDD-9CB65404EEF0.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/EC739076-5819-E411-9B94-001E673C84B9.root',
       '/store/mc/Spring14miniaod/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/F49687A8-5919-E411-AD79-001E6739CEB1.root' ] );


secFiles.extend( [
               ] )

