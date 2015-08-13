import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/2E3CB549-2A09-E411-9EC0-20CF30561706.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/32943B58-2A09-E411-B75F-52540097313B.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/40012D14-2A09-E411-835A-0025907B50D6.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/44167F5D-2A09-E411-B2C3-0025907B4FE4.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/6E8FD94B-2A09-E411-AD66-20CF3027A57B.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/841C3C14-2A09-E411-BD4B-20CF305B04DA.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/92C04017-2A09-E411-BCE6-002590D0B0C0.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A4B01A40-2A09-E411-8656-002590D0B000.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E052DF42-2A09-E411-A1EB-0025907B4EF0.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F294493E-2A09-E411-ABEE-20CF30561701.root',
       '/store/mc/Spring14miniaod/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/FEED0446-2A09-E411-B3AA-002590760990.root' ] );


secFiles.extend( [
               ] )
