import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0018DF4B-A509-E411-8FE8-0024E85A4BF4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0E05E175-A609-E411-B3C3-0024E85A4BF4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/106A86AC-A409-E411-925C-0024E85A3F71.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/120CB6FA-B409-E411-814F-00221981AF26.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/16310E0F-B109-E411-9CF8-0024E856F86E.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/206C3C56-A609-E411-9271-68B599B9E940.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/22DE772D-B109-E411-B2C3-00221980E840.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/26B4AAB0-9909-E411-8DD6-002219826F48.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/28190DBD-B409-E411-B488-0024E850DF94.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/2E8E4E49-A809-E411-8D40-0024E850DF94.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/301F2CE1-2709-E411-984A-0024E85A405A.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3A7B3E8D-AD09-E411-AAF8-002219826BCD.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3A99235E-AC09-E411-B9B6-0024E85A4092.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3ACF96DF-A909-E411-83DF-002219818E18.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3CC57176-B809-E411-B19B-00221981BAA3.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/3CD58E51-AB09-E411-BE3B-0024E850DF90.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4089C0CD-A909-E411-B848-0024E850DF90.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/4C871FA7-A309-E411-A5E4-0024E85A4BF4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/5808294C-A909-E411-BEC7-002219818E18.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/64C2A7F6-B409-E411-A049-0024E850DF94.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/741693C4-B509-E411-B1B5-14FEB5FB6872.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7C2CCA4D-A509-E411-A470-0024E85A3EE4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/7E5D509F-B809-E411-B944-00221981BAA3.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8076DE76-A609-E411-AB0F-0024E85A3EE4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/82B1FCDB-9F09-E411-9B42-0024E85A405A.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/84B3C5A5-9F09-E411-A258-0024E85A4FB7.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8AF7F8C2-B409-E411-B575-002219826F48.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8AFC48C0-A309-E411-9F3C-0024E85A3EE4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8CA4015C-A509-E411-A49F-00221981B44C.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/8EB67B7C-9C09-E411-8CC2-0024E850DF90.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/962E1791-9909-E411-8F96-0024E85A4BFC.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/9C95792C-B109-E411-BDE2-00221981AF26.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A0701A59-A509-E411-BBB0-0024E85A4BEC.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/AC0682B1-AC09-E411-9CC0-0024E85A4C00.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/ACA8A329-9E09-E411-A307-0024E850DF81.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B4E8DFEC-B509-E411-9FC5-0024E850DF94.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B8CB9AC8-B509-E411-9207-0024E850DF94.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BA3B5C41-BE09-E411-B702-0024E850DF90.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/BA479E2D-B109-E411-AD61-00221980E840.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C0F39852-A809-E411-A0E5-0024E856F86E.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D6049E4D-A909-E411-9E66-0024E85A4092.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/D8F53E5D-CE09-E411-812F-00221981BAAB.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/DCD7FE57-9C09-E411-9BC3-0024E85A08F8.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E04361AB-A309-E411-AF4C-0022198273D4.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E0794AB0-CA09-E411-937E-002219826F44.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E0D259D1-AA09-E411-AE8C-0024E85A4092.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E212CCFF-BD09-E411-BBBC-0024E85A4BEC.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/EA4EDB24-A009-E411-961D-00221981AF36.root',
       '/store/mc/Spring14miniaod/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F4C1D713-C909-E411-8400-68B599B9D9B8.root' ] );


secFiles.extend( [
               ] )

