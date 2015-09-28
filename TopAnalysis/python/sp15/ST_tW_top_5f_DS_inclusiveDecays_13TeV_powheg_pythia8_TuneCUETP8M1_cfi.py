import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/02D835D8-FC01-E511-A105-525400929F2C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/22894C59-FD01-E511-AAA1-E0CB4E553643.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/24F64413-EB01-E511-98B9-0025905964A6.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/2665ED5A-FD01-E511-A685-00259073E3F0.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/46E51013-EB01-E511-9136-0026189438C9.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/72A37118-EB01-E511-9FEA-0025905A609E.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/9C5D08DC-FC01-E511-891F-0025907B4EEA.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/A06764A2-EA01-E511-AE7D-0025B3E0652A.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/CE87C51A-EB01-E511-9C92-0025905A60A6.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/D218475A-EA01-E511-89A6-002590200B38.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/EA1F2233-EA01-E511-A746-002590A88736.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/EE472C53-FD01-E511-B8D6-20CF3027A5C4.root' ] );


secFiles.extend( [
               ] )


