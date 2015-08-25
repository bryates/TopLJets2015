import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/22765F90-46FF-E411-AB4A-00259073E438.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/2CAD31D2-46FF-E411-841E-00074305CC86.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/4073F5F7-46FF-E411-93F8-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/48DA79C8-46FF-E411-BFCD-0025905A60B8.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/6E11FBCB-46FF-E411-B5B5-002590574604.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/0A388965-48FF-E411-BC9F-00259077501E.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/3267E471-48FF-E411-9C20-00074305CF01.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/345A4067-48FF-E411-A53D-00074305CB9B.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/704C447B-48FF-E411-9A4F-00074305CCA1.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/803D5E65-48FF-E411-895C-00259073E410.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/A4D65E6B-48FF-E411-8E7C-0025905A6092.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/B8DE447D-48FF-E411-A100-00074305CB90.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/C09E4465-48FF-E411-AB24-00259073E466.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/CEB61569-48FF-E411-87FE-00074305CB49.root',
       '/store/mc/RunIISpring15DR74/ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/D0922670-48FF-E411-8601-0025905B855C.root' ] );


secFiles.extend( [
               ] )


