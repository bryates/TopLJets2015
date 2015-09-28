import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/00839933-6202-E511-A035-0025904A8EC4.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/044EF786-4602-E511-B010-0025905A60DE.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/0CCE5DDC-6102-E511-AF94-0CC47A4DEDCA.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/167711A0-4602-E511-8689-0026189438CB.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/4E1830A8-4602-E511-97EB-00259059391E.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/687746D6-FC01-E511-86DC-0025907B50D2.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/6C95AA2B-6202-E511-A244-0025905B857C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/70CAFAE0-1902-E511-BAA4-0025907254C8.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/84884A45-FD01-E511-BE90-002590D0AFF2.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/8AE757FF-6102-E511-85CB-00074305CE1D.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/8CF5112B-6202-E511-88AE-0025905A48C0.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/9012B3B9-FD01-E511-A949-E0CB4E55367D.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/A08BE84F-F401-E511-8367-A0040220FE80.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/A60E8414-4702-E511-BBBA-0025905A60B8.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/A66BD52C-FE01-E511-81C4-002590574A44.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/AC01DEB5-FD01-E511-BC29-20CF3027A629.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/AECA4E29-FE01-E511-A405-0CC47A4DEE70.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/B0DB62A5-4602-E511-B9A5-0025905A60B6.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/BE99DC11-DB01-E511-B8E1-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/C08FD72D-FE01-E511-9FEE-00259073E4B6.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/D01E1222-DB01-E511-A271-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/F06ACE8B-FD01-E511-B304-002590D0AFD2.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/F2335681-4602-E511-B06C-0025905A60B0.root' ] );


secFiles.extend( [
               ] )


