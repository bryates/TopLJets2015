import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0A2744F9-FA05-E411-BD0C-00259073E36C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0E936434-FD05-E411-81BF-F4CE46B27A1A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/32E07232-FD05-E411-897C-00259073E522.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3CE2B535-FB05-E411-919A-20CF307C98DC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/48093276-FC05-E411-9EEE-001F296564C6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/50B66FF3-FA05-E411-A937-001F296564C6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/544B2DF7-FA05-E411-B91F-001F2965F296.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/54DB2FF7-FE05-E411-824B-00259073E522.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/56D1BC32-FD05-E411-A512-20CF3027A5EB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5AD70432-FC05-E411-906C-20CF3027A5CD.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5C4FBFF4-FA05-E411-9767-00259073E36C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5CF748F8-FC05-E411-814B-20CF3027A5A2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7806E24D-FC05-E411-8922-001F2965F296.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7C16B231-FD05-E411-8E00-20CF3027A5EB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/802452C1-FC05-E411-A969-00221983E092.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8217E3BD-FC05-E411-B8C2-0025907277CE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8676BEF4-FA05-E411-B26A-00259073E36C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8C1741F3-FA05-E411-B5B5-20CF3027A582.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8C915AB8-FC05-E411-9EAF-F4CE46B27A1A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AA0FCBB0-FC05-E411-898D-00259073E36C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B49383BA-FC05-E411-9914-F4CE46B27A1A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B6DAEFDD-FB05-E411-9851-20CF3027A5CD.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C6F5C44F-FD05-E411-B86F-D48564592B02.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C83B6B6C-FC05-E411-BAFD-D48564599CAA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/CEF64C64-FD05-E411-A799-001F2965648A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D6C305FC-FA05-E411-9AF5-00259073E522.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/DE0FC6A4-FC05-E411-A2F9-00259073E36C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E2D5AD33-FD05-E411-868A-D48564594F36.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E63BCC43-FB05-E411-834E-D48564599CEE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/EAD01F32-FD05-E411-91E4-20CF3027A5F4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F0A18D25-FC05-E411-8DFC-20CF3027A582.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F0B8E6B6-FA05-E411-9DAE-20CF3027A5CD.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F23A21C3-FD05-E411-9E29-A4BADB3D00FF.root' ] );


secFiles.extend( [
               ] )

