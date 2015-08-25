import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/628B3F1A-28FF-E411-BEDF-20CF3019DF09.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/6CF47F80-27FF-E411-B98E-002590A831B4.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/B24D4580-27FF-E411-A7AE-002590200964.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/489BD5F5-10FF-E411-B8E8-001E67396E3C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/F6D9ABFB-10FF-E411-9F91-E0CB4E29C4D7.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/FE2AD6FA-10FF-E411-B52F-00259073E4A0.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/06BCC94F-E9FE-E411-8ECD-485B39800B67.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/3E4D3796-E9FE-E411-B343-0025B3E0656C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/4C670393-E9FE-E411-BDF0-002590A371AE.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/76B37A50-E9FE-E411-8FDF-00259077501E.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/DA9E274F-E9FE-E411-AE7D-00259073E514.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/324B0854-04FF-E411-A943-E0CB4E29C4DB.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/4C42C0FD-04FF-E411-9847-0025905B858C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/5631B8FD-04FF-E411-9133-0025905B858C.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/C01FB4AE-04FF-E411-AF74-001E67396707.root',
       '/store/mc/RunIISpring15DR74/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/DC1ECC58-04FF-E411-9575-BCAEC51B8F58.root' ] );


secFiles.extend( [
               ] )


