import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/00F3C52D-FB05-E411-B149-002590AC4BF6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/08858AAC-FA05-E411-97CA-0025904B5FB8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1C7E39AF-FA05-E411-8159-002590494C92.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/36885A14-FB05-E411-BF5C-002590AC4FC8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4422FED7-FA05-E411-8684-002590494E34.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/48C9E982-FA05-E411-A1C4-002481E0D6EE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4E71ADAC-FA05-E411-9BD8-003048F30422.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/543EFDCF-FA05-E411-BDB6-002590494E36.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/64C53807-FB05-E411-A89F-00266CF32F14.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/68F16C25-FB05-E411-BB7F-003048CEFFE4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8EC2AC85-FA05-E411-A392-0025904B12A8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A6BE1922-FB05-E411-AE2B-00266CFFA0B0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AA29663C-A305-E411-A49C-002590494E38.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AE747724-FB05-E411-93C0-D8D385FF6C5E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B061CFA1-FA05-E411-8279-002590AC5482.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B070E92C-FB05-E411-873E-002590AC4C6C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B60156C0-FA05-E411-A047-002481E94BFE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B6BFC528-FB05-E411-9DBA-0025904B0FB6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BA2F03B4-FA05-E411-A11B-0025901D493E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D034FC26-FB05-E411-BC7C-0025904B1446.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/DE80CFAE-FA05-E411-8025-00266CFFA048.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E49E4014-FB05-E411-960F-002590AC4C74.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E6F5D348-A305-E411-9B1F-0025901D4090.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/EE38EF11-FB05-E411-9723-002590AC4C6E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/EE8BEA09-FB05-E411-A02E-002481E0DA96.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F63E61B0-FA05-E411-A94A-0025901D4C94.root',
       '/store/mc/Spring14miniaod/QCD_Pt-170to300_MuEnrichedPt5_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/FCE8D58B-FA05-E411-A6EC-0025907FD3CC.root' ] );


secFiles.extend( [
               ] )

