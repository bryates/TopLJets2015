#!/usr/bin/env python

import ROOT
from runWeightCreationLoop import TUNES,MAXEVENTS

BRs=[
    (511, 'B^{0}',        0.23845, 0.1043, 0.1033, 0.0028),
    (521, 'B^{+}',        0.25579, 0.1129, 0.1099, 0.0028),
    (531, 'B^{0}_{s}',    0.21920, 0.0930, 0.0960, 0.008),
    (5122, '#Lambda_{b}', 0.17870, 0.0770, 0.104, 0.022),
    ('inc', 'inc', 1, 1, 1, 1)
    ]

"""
Interpolate extremes and then derive the weights based on a 2nd order spline for the remaining nodes
"""
def smoothWeights(gr):

    #interpolate for low xb
    gr.Fit('pol3','QR+','',0,0.7)
    lowxb=gr.GetFunction('pol3')
    #gr.Fit('pol9','QR+','',0,0.7)
    #lowxb=gr.GetFunction('pol9')

    #flatten tail for xb>1
    gr.Fit('pol0','QR+','',1.03,2)
    highxb=gr.GetFunction('pol0')
    if 'Dz413' in gr.GetName() or 'Dzrand413' in gr.GetName(): 
	#gr.Draw('ap')
        gr.Fit('pol4','QR','',0,0.8)
	lowxb=gr.GetFunction('pol4')
	gr.Fit('pol0', 'QR+', '', 0.6, 1)
	highxb=gr.GetFunction('pol0')

    smoothExtremesGr=ROOT.TGraph()
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,gr.GetN()):
        gr.GetPoint(i,x,y)
        if 'Dz413' in gr.GetName() or 'Dzrand413' in gr.GetName(): 
            if x<0.8:
                y_fit = 1.05639 + 1.57986 * pow(x,1) + -15.3383 * pow(x,2) + 34.9746 * pow(x,3) + -23.3951 * pow(x,4)
	        y_fit = 1.18739 + 1.31977 * pow(x,1) + -19.3186 * pow(x,2) + 45.6512 * pow(x,3) + -30.1078 * pow(x,4)
            else:
                y_fit =  1.01233
	        y_fit = 1.02571
            #if x<0.1:
            #    y_fit = 1.05435 + 5.35532 * pow(x,1) + -43.1835 * pow(x,2)
            #elif x<0.8:
            #    y_fit = 1.72476 + -8.55134 * pow(x,1) + 32.806 * pow(x,2) + -51.8687 * pow(x,3) + 28.7536 * pow(x,4) + 1.77362 * pow(x,5)
            #else:
            #    y_fit = 2.24065
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzuinc' in gr.GetName():
            if x<0.6:
                y_fit = 1.00305 + 0.193647 * pow(x,1) + -0.919419 * pow(x,2) + 0.908652 * pow(x,3)
            if x<1.0:
	        y_fit = 0.899912 + 0.36203 * pow(x,1) + -0.383503 * pow(x,2)
            else:
	        y_fit = 0.904609
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzdinc' in gr.GetName():
            if x<0.6:
		y_fit = 1.00114 + -0.195864 * pow(x,1) + 0.861775 * pow(x,2) + -0.838063 * pow(x,3)
            elif x<1.0:
	        y_fit = 1.07757 + -0.271801 * pow(x,1) + 0.285673 * pow(x,2)
            else:
                y_fit = 1.05969
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzuchargedinc' in gr.GetName():
            if x<0.6:
	        y_fit = 1.00532 + 0.0351235 * pow(x,1) + -0.0425894 * pow(x,2) + -0.0579174 * pow(x,3)
            elif x<1.0:
		y_fit = 1.02612 + -0.0522146 * pow(x,1) + 0.0148524 * pow(x,2)
            else:
	        y_fit = 0.989387
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzdchargedinc' in gr.GetName():
            if x<0.6:
                y_fit = 1.00702 + -0.142328 * pow(x,1) + 0.384509 * pow(x,2) + -0.283077 * pow(x,3)
            elif x<1.0:
	        y_fit = 0.983544 + 0.0288486 * pow(x,1) + -0.00473737 * pow(x,2)
            else:
	        y_fit = 1.00764
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzdown' in gr.GetName():
            if x<0.6:
	        y_fit = 1.00007 + -0.147383 * pow(x,1) + 0.351221 * pow(x,2) + -0.17774 * pow(x,3)
            elif x<1.0:
	        y_fit = 0.957158 + 0.0829607 * pow(x,1) + -0.0224346 * pow(x,2)
            else:
	        y_fit = 1.02168
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue
        elif 'Dzup' in gr.GetName():
            if x<0.6:
		y_fit = 0.992566 + 0.5172 * pow(x,1) + -1.35814 * pow(x,2) + 0.867965 * pow(x,3)
            elif x<1.0:
		y_fit = 1.10619 + -0.192541 * pow(x,1) + 0.0353941 * pow(x,2)
            else:
	        y_fit = 0.950417
	    smoothExtremesGr.SetPoint(i,x,y_fit)
            continue


        if x<0.55:
            smoothExtremesGr.SetPoint(i,x,lowxb.Eval(x))
        elif x>1.03:
            smoothExtremesGr.SetPoint(i,x,highxb.Eval(x))
        else:
            smoothExtremesGr.SetPoint(i,x,y)

    #smooth the weights
    tSpline=ROOT.TMVA.TSpline2("spline",smoothExtremesGr)
    smoothGr=ROOT.TGraph()
    for i in xrange(0,1000):
        x=2.0*i/1000
        smoothGr.SetPoint(i,x,tSpline.Eval(x))
    return smoothGr


def main():

    outf='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/bfragweights.root'
    
    #derive the weights
    xb={}
    for tag,_,_ in TUNES:
        if ('Dss' in tag or 'Dz' in tag) and 'charged' not in tag:
            fName='LJets2015/2016/bfrag/xb_%s_numEvent%d.root'%('fit',MAXEVENTS) if MAXEVENTS>0 else 'LJets2015/2016/bfrag/xb_%s.root'%'fit'
        else:
            fName='LJets2015/2016/bfrag/xb_%s_numEvent%d.root'%(tag,MAXEVENTS) if MAXEVENTS>0 else 'LJets2015/2016/bfrag/xb_%s.root'%tag
        #fIn=ROOT.TFile.Open(fName)
        fIn=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_%s.root'%tag)
        if 'Dz' in tag:
	    fIn=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fit.root')
        #xb[tag]=fIn.Get('bfragAnalysis/xb_semilepinc').Clone(tag)
        if 'rand' in tag:
            print(tag)
            fIn2=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fit.root')
	    xb[tag.replace('rand', '')]=fIn2.Get('bfragAnalysis/xb_%s' % tag.replace('rand', '')).Clone()
	    xb[tag]=fIn2.Get('bfragAnalysis/xb_%s' % tag).Clone()
            print(xb[tag].Integral())
        elif any([n in tag for n in ['DssDz','431','413']]) and any([n in tag for n in ['u','d']]):
            xb[tag]=fIn.Get('bfragAnalysis/xb_inc').Clone(tag)
            tmp=fIn.Get('bfragAnalysis/xb_%s' % tag).Clone(tag+'tmp')
            fIn2 = []
            if 'Dzu' in tag:
                fIn2=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fitup.root')
            elif 'Dzd' in tag:
                fIn2=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fitdown.root')
            else:
                pass
            print('xb_'+tag)
            tmpD=fIn2.Get('bfragAnalysis/xb_%s' % tag).Clone(tag+'tmpD')
            xb[tag].Add(tmp, -1) # Remove e.g. DssDz
            xb[tag].Add(tmpD, 1) # Add in e.g. DssDz shifted up
            fIn2.Close()
        #elif 'DssDz' == tag:
        #    xb[tag]=fIn.Get('bfragAnalysis/xb_%s' % tag).Clone()#' % tag).Clone()
        elif 'Dss' in tag and 'noDss' not in tag:
            xb[tag]=fIn.Get('bfragAnalysis/xb_inc').Clone(tag)
            if any(n == tag[-1] for n in ['u','d']):
                tmp=fIn.Get('bfragAnalysis/xb_%s' % tag[:-1]).Clone()#' % tag).Clone()
                xb['DssDz'] = fIn.Get('bfragAnalysis/xb_%s' % tag[:-1]).Clone()#' % tag).Clone()
                if tag[-1] == 'u':
                    fIn2=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fitup.root')
                else:
                    fIn2=ROOT.TFile.Open('LJets2015/2016/bfrag/xb_fitdown.root')
                tmpD=fIn2.Get('bfragAnalysis/xb_%s' % tag[:-1]).Clone()#' % tag).Clone()
                xb[tag].Add(tmp, -1)
                xb[tag].Add(tmpD, 1)
            elif 'Ten' in tag:
                pass
                tmp=fIn.Get('bfragAnalysis/xb_%s' % 'DssDz').Clone()#' % tag).Clone()
                xb[tag].Add(tmp, -1)
                xb[tag + 'TenUp'] = xb[tag].Clone()
                xb[tag + 'TenDown'] = xb[tag].Clone()
                tmp.Scale(1.1) # Vary Dss by 10%
                xb[tag + 'TenUp'].Add(tmp) # Vary Dss by 10%
                tmp.Scale(1/1.1) # Undo vary Dss by 10%
                tmp.Scale(0.9) # Vary Dss by -10%
                xb[tag + 'TenDown'].Add(tmp) # Vary Dss by 10%
            else:
                pass
        #elif any(n in tag for n in ['BpmDz', 'BsDz', 'LbDz']):
        #    d = {'Bpm' : 521, 'Bs' : 531, 'Lb' : 5122}
        #    xb[tag]=fIn.Get('bfragAnalysis/xb_Dzu%s' % d[tag[:-2]]).Clone(tag)
        #    xb['B0Dz']=fIn.Get('bfragAnalysis/xb_Dzu511').Clone('B0Dz')
        else:
            xb[tag]=fIn.Get('bfragAnalysis/xb_inc').Clone(tag)
        #xb[tag].Rebin(2);
        xb[tag].SetDirectory(0)
        fIn.Close()

    #save to file
    fOut=ROOT.TFile.Open(outf,'RECREATE')
    #for tag in ['up','uup','uuup','central','ccentral','cccentral','down','ddown','dddown']:
    #for tag in [ ('sup','BL',0.865), ('scentral','BL',0.854), ('sdown','BL',0.845) ]
    for tag,frag,param in TUNES:
    #for tag in ['up','central','down']:
    #for tag in ['up','central','down','Peterson']:
        if tag == 'cuetp8m2t4': continue
        if 'Bhadron' in tag: continue
        print(tag)
        xb[tag].Scale(1./xb[tag].Integral())
        xb['cuetp8m2t4'].Scale(1./xb['cuetp8m2t4'].Integral())
        if 'rand' in tag:
            xb[tag.replace('rand', '')].Scale(1./xb[tag.replace('rand', '')].Integral())
            xb[tag].Divide(xb[tag], xb[tag.replace('rand', '')], 1./xb[tag].Integral(), 1./xb[tag.replace('rand', '')].Integral())
        else: xb[tag].Divide(xb['cuetp8m2t4'])
        #xb[tag].Smooth()
        gr=ROOT.TGraphErrors(xb[tag])
        gr.SetName(tag)
        gr.SetMarkerStyle(20)

        sgr=smoothWeights(gr)
        sgr.SetName(tag+'Frag')
        sgr.SetLineColor(ROOT.kBlue)
        sgr.Write()

        #gr.Draw('ap')
        #sgr.Draw('l')
        #raw_input()

    #semi-leptonic BRs
    semilepbrUp=ROOT.TGraph()
    semilepbrUp.SetName("semilepbrUp")
    semilepbrDown=ROOT.TGraph()
    semilepbrDown.SetName("semilepbrDown")

    for entry in BRs:

        i=semilepbrUp.GetN()
        pid,_,py8inc,py8exc,pdg,pdgUnc=entry
        if pid == 'inc': continue

        brUp   = py8inc*(1+ROOT.TMath.Max((pdg+pdgUnc)-py8exc,0.)/py8exc)
        semilepbrUp.SetPoint(i,     pid*(-1), (1-brUp)/(1-py8inc))
        semilepbrUp.SetPoint(i+1,   pid,      brUp/py8inc)

        brDown = py8inc*(1-ROOT.TMath.Max(py8exc-(pdg-pdgUnc),0.)/py8exc)
        semilepbrDown.SetPoint(i,   pid*(-1), (1-brDown)/(1-py8inc))
        semilepbrDown.SetPoint(i+1, pid,      brDown/py8inc)

    semilepbrUp.Write()
    semilepbrDown.Write()

    #for entry in BRs:
        #for unc in ['Dzb','Dzbu','Dzbd']:#,'Dz','Dzu','Dzd']:
            #pid=entry[0]
            #if pid == 'inc' and unc[-1] == 'u': continue
            #if pid == 'inc' and unc[-1] == 'd': continue
            #if pid != 'inc' and unc[-1] == 'z': continue
            #if pid != 'inc' and unc[-1] == 'b': continue
#
            #name='bfragAnalysis/xb_semilep{}{}'.format(unc,pid)
            #fIn=ROOT.TFile.Open(fName)
            #tag=unc+str(pid)
            #xb[tag]=fIn.Get(name).Clone(unc+str(pid))
            #if xb[tag].Integral():
                #xb[tag].Scale(1./xb[tag].Integral())
            #xb[tag].Divide(xb['cuetp8m2t4'])
            #xb[tag].SetDirectory(0)
            #fIn.Close()

    #for entry in BRs:
        #if entry[0] == 'inc': continue
        #for unc in ['Dzb']:#,'Dz']:
            #pid=entry[0]
            #tag=unc+str(pid)
            #tagu=unc+'u'+str(pid)
            #tagd=unc+'d'+str(pid)
            #x = xb[tagd].Clone(tag+'over'+tagd)
            ##x.Divide(x,xb[tag],1,1,'B')
            #if x.Integral(): x.Scale(1./x.Integral())
            #x.Divide(xb['cuetp8m2t4'])
            #fOut.cd()
            #x.Write()
            #gr=ROOT.TGraphErrors(x)
            #gr.SetMarkerStyle(20)
            #sgr=smoothWeights(gr)
            #sgr.SetName('semilep'+unc+str(pid)+'Frag')
            #sgr.SetLineColor(ROOT.kBlue)
            #fOut.cd()
            #sgr.Write()
#
    #for unc in ['Dzb']:#,'Dz']:
        #pid=entry[0]
        ##xb_semilepDzbinc
        #tag='xb_semilep'+unc+'inc'
        #tag=unc+'inc'
        #tagu=unc+'u'+str(pid)
        #tagd=unc+'d'+str(pid)
        #name='bfragAnalysis/xb_semilep{}{}'.format(unc,pid)
        #print name
        ##x=fIn.Get(name).Clone(unc+str(pid))
        #x = xb[tag].Clone(tag+'overnom')
        ##x.Divide(x,xb[tag],1,1,'B')
        ##x.Scale(1./x.Integral())
        ##x.Divide(xb['cuetp8m2t4'])
        ##x.Divide(x,xb['cuetp8m2t4'],1,1,'B')
        #fOut.cd()
        #x.Write()
        #gr=ROOT.TGraphErrors(x)
        #gr.SetMarkerStyle(20)
        #sgr=smoothWeights(gr)
        #sgr.SetName('semilep'+unc+'uncFrag')
        #print sgr.GetName()
        #sgr.SetLineColor(ROOT.kBlue)
        #fOut.cd()
        #sgr.Write()

    fOut.Close()
    print 'Fragmentation been saved to',outf

if __name__ == "__main__":
    main()
