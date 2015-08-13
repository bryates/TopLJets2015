import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/107FE127-A806-E411-8857-0025904956BE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/16E0941F-AA06-E411-8757-00259074AE34.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/2ACC80F0-A606-E411-AB34-20CF3027A668.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4CFB0B7A-C306-E411-A1C5-20CF3027A686.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/50AD5CFA-A006-E411-A448-20CF3027A668.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/78C42EDB-B806-E411-8202-20CF305B066B.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/80438322-AA06-E411-BCCF-00259074AE34.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A0FBC8CD-A906-E411-A2E6-20CF307C98F1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A8C37D37-CF06-E411-A181-00259074AEA2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AE8CBCAA-BC06-E411-808B-20CF305B066B.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C00B1303-EB06-E411-B440-20CF30725213.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/CE18C9CD-A906-E411-A49C-20CF307C98F1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D80E5C51-C006-E411-B6E6-20CF305B066B.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E2D26620-AA06-E411-B5C8-002590D4FBDA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E8E5E62C-AA06-E411-9BDC-20CF3027A668.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F4BA7620-AA06-E411-AABC-002590D4FBDA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-600to800_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F8C763BD-AA06-E411-B835-002590D4FBDA.root' ] );


secFiles.extend( [
               ] )

