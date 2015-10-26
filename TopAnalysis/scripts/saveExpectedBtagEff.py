import ROOT
from array import array
from TopLJets2015.TopAnalysis.storeTools import *

csvWP='j_csv>0.890'
inputDir='/store/cmst3/user/psilva/LJets2015/7bae03e/MC13TeV_TTJets'

#open a ttbar file
input_list=getEOSlslist(directory=inputDir)
data=ROOT.TChain('analysis/data')
#for i in xrange(0,len(input_list)):
for i in xrange(0,1):
    data.Add(input_list[i])

print 'Projecting tagging efficiency from ',data.GetEntries(),' events'
print 'Working point is',csvWP


#prepare histograms
effgrs=[]
ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV/c];',len(ptBins)-1,array('d',ptBins))
preTagH.Sumw2()
tagH=preTagH.Clone('tagH')

#count number of tagged jets
for flav,cond in [('b',"abs(j_hadflav)==5"),
                  ('c',"abs(j_hadflav)==4"),
                  ('udsg','abs(j_hadflav)!=5 && abs(j_hadflav)!=4')]:
    print '\t computing for',flav,cond
    preTagH.Reset('ICE')
    tagH.Reset('ICE')
    data.Draw("j_pt >> preTagH",cond,'goff')
    data.Draw('j_pt >> tagH',cond + ' && ' + csvWP,'goff')
    effgrs.append(ROOT.TGraphAsymmErrors())
    effgrs[-1].SetName(flav)
    effgrs[-1].Divide(tagH,preTagH)

#save efficiency graphs
fOut=ROOT.TFile.Open('data/expTageff.root','RECREATE')
for gr in effgrs:
    gr.Write()
fOut.Close()
