import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0E6D2ECF-FA05-E411-9FF8-00E0814189CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/46ED6304-FB05-E411-B02E-00E08133C86B.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4C02A0C6-FA05-E411-B326-002618943C0A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/52AAB6B0-FB05-E411-BAE0-002618943C0D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7A4F75DA-FA05-E411-AF20-BCAEC54E98B4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7C7C5387-FA05-E411-BAF3-002618943C0D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7C979F70-FA05-E411-82C8-002618943C0A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/9C653F69-FB05-E411-BCB2-00E08132728F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BC4993D5-FA05-E411-BE0F-F46D0450CEA0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C46A9069-FA05-E411-83CD-3085A9262DB8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C61FF9C6-FA05-E411-BDDC-BCAEC51FDEED.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C64FB9C2-FA05-E411-B315-F46D0450CE46.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C84BE4C6-FA05-E411-9847-3085A9262DB8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/CAA89CC8-FA05-E411-B83A-002618943C31.root',
       '/store/mc/Spring14miniaod/QCD_Pt-800to1000_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E4C950C5-FA05-E411-9146-3085A9262DA0.root' ] );


secFiles.extend( [
               ] )

