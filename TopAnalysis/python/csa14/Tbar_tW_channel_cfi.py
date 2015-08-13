import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/0C3A0A4B-E318-E411-B033-002590D0AFFC.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/105452AC-E418-E411-995E-E0CB4E29C4F9.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/2E095685-E418-E411-A703-20CF300E9ECF.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/2E235F31-E318-E411-A6F4-E0CB4E29C4F9.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/4E70ED20-E418-E411-BEFB-0025907B4F28.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/522C32EC-E118-E411-BF7F-00259073E4C2.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/5299C196-E518-E411-AE8A-E0CB4EA0A92E.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/56109483-E218-E411-BBAB-E0CB4EA0A92E.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/56463A6A-E318-E411-945D-00259073E4C2.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/607D5313-E518-E411-AC01-00259073E44C.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/7025070D-E418-E411-8C90-20CF3027A5D5.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/7E9782C3-E118-E411-B73D-E0CB4E29C4F9.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/8608261C-E518-E411-8AA8-00259073E390.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/8CA5E813-E218-E411-A980-002590D0AFFC.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/8CDAFA0B-E418-E411-9A37-E0CB4EA0A92E.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/A8A9ECD8-E418-E411-919F-00259073E4C2.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/B462E69D-E318-E411-BE1C-00259073E37C.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/D8F9780D-E518-E411-AB04-00259073E37C.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/DADCFFB5-E418-E411-A734-002590D0AFFC.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/EA4F259E-E518-E411-A428-0025907B4F28.root',
       '/store/mc/Spring14miniaod/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/30000/F896A64A-E518-E411-99DE-20CF3027A5D5.root' ] );


secFiles.extend( [
               ] )

