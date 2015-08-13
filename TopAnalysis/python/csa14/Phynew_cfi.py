import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1C772C5B-DD71-E411-B64E-0025904B1370.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2010CF0F-CB71-E411-9331-002481E0D678.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/20A8A674-DD71-E411-A04B-0025901D4B22.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5042A9C5-C971-E411-885B-0025901D42C0.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64B14A09-CB71-E411-98A3-0025907FD24A.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2C09A49-CB71-E411-9ACE-002590494C92.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F697A0CD-C971-E411-9404-003048D3CB3C.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C2923647-C671-E411-86E0-003048D438FE.root',
       '/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E2C2C43B-C671-E411-A551-002590AC4B76.root' ] );


secFiles.extend( [
               ] )



