import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle

"""
test if variation is significant enough i.e. if sum_{bins} |var-nom| > tolerance
"""
def acceptVariationForDataCard(nomH,upH,downH,tol=1e-2):
    diffUp,diffDown=0,0
    for xbin in xrange(1,nomH.GetNbinsX()):
        val,valUp,valDown=nomH.GetBinContent(xbin),upH.GetBinContent(xbin),downH.GetBinContent(xbin)
        diffUp+=ROOT.TMath.Abs(valUp-val)
        diffDown+=ROOT.TMath.Abs(valDown-val)
    accept = True if (diffUp>tol or diffDown>tol) else False
    return accept


"""
get distributions from file
"""
def getDistsFrom(directory):
    obs=None
    exp={}
    dirName=directory.GetName()
    for key in directory.GetListOfKeys():
        obj=directory.Get(key.GetName())
        if not obj.InheritsFrom('TH1') : continue
        if obj.GetName()==dirName : 
            obs=obj.Clone('data_obs')
            obs.SetDirectory(0)
        else : 
            newName=obj.GetName().split(dirName+'_')[-1]
            for token in ['+','-','*',' ','#','{','(',')','}','@']:
                newName=newName.replace(token,'')
            exp[newName]=obj.Clone(newName)
            exp[newName].SetDirectory(0)
            if exp[newName].InheritsFrom('TH2'):
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    for ybin in xrange(1,exp[newName].GetNbinsY()+1):
                        binContent=exp[newName].GetBinContent(xbin,ybin)
                        if binContent>0: continue
                        newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                        exp[newName].SetBinContent(xbin,ybin,newBinContent)
                        exp[newName].SetBinError(xbin,ybin,newBinContent)
            else:
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    binContent=exp[newName].GetBinContent(xbin)
                    if binContent>0: continue
                    newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                    exp[newName].SetBinContent(xbin,newBinContent)
                    exp[newName].SetBinError(xbin,newBinContent)
    return obs,exp

"""
save distributions to file
"""
def saveToShapesFile(outFile,shapeColl,directory=''):
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    if len(directory)==0:
        fOut.cd()     
    else:
        outDir=fOut.mkdir(directory)
        outDir.cd()
    for key in shapeColl:
        #remove bin labels  
        shapeColl[key].GetXaxis().Clear()

        #convert to TH1D (projections are TH1D)
        if not shapeColl[key].InheritsFrom('TH1D') :
            h=ROOT.TH1D()
            shapeColl[key].Copy(h)
            shapeColl[key]=h

        shapeColl[key].Write(key,ROOT.TObject.kOverwrite)
    fOut.Close()


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',          dest='input',       help='input plotter',                                      default=None,          type='string')
    parser.add_option(      '--systInput',      dest='systInput',   help='input plotter for systs from alternative samples',   default=None,          type='string')
    parser.add_option('-d', '--dist',           dest='dist',        help='distribution',                                       default='njetsnbtags', type='string')
    parser.add_option('-s', '--signal',         dest='signal',      help='signal (csv)',                                       default='tbart',       type='string')
    parser.add_option('-c', '--cat',            dest='cat',         help='categories (csv)',                                   default='1j,2j,3j,4j', type='string')
    parser.add_option('-q', '--qcd',            dest='qcd',         help='qcd normalization file',                             default=None,          type='string')
    parser.add_option('-o', '--output',         dest='output',      help='output directory',                                   default='datacards',   type='string')
    (opt, args) = parser.parse_args()

    signalList=opt.signal.split(',')
    catList=opt.cat.split(',')

    #read qcd normalization
    qcdNormUnc=None
    if opt.qcd:
        cache=open(opt.qcd,'r')
        pickle.load(cache)
        qcdNormUnc=pickle.load(cache)

    #prepare output directory 
    os.system('mkdir -p %s'%opt.output)

    #get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input)
    systfIn=None
    if opt.systInput:
        systfIn=ROOT.TFile.Open(opt.systInput)

    #loop over categories
    for cat in catList:

        print 'Iniating datacard for ',cat

        #nomimal expectations
        obs,exp=getDistsFrom(directory=fIn.Get('%s_%s'%(opt.dist,cat)))
    
        #prepare output ROOT file
        outFile='%s/shapes_%s.root'%(opt.output,cat)
        fOut=ROOT.TFile.Open(outFile,'RECREATE')
        fOut.Close()

        #start the datacard
        datacard=open('%s/datacard_%s.dat'%(opt.output,cat),'w')
        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis category %s\n' % (getpass.getuser(),
                                                                                          commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                                                                                          cat) )

        datacard.write('#\n')
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')
        datacard.write('-'*50+'\n')
        datacard.write('shapes *        * shapes_%s.root nom/$PROCESS $SYSTEMATIC/$PROCESS\n' % cat)
        datacard.write('-'*50+'\n')
        datacard.write('bin 1\n')
        datacard.write('observation %3.1f\n' % obs.Integral())
        datacard.write('-'*50+'\n')

        #expectations
        datacard.write('\t\t\t %15s'%'bin')
        for i in xrange(0,len(exp)): datacard.write('%15s'%'1')
        datacard.write('\n')
        datacard.write('\t\t\t %15s'%'process')
        for sig in signalList: datacard.write('%15s'%sig)
        for proc in exp: 
            if proc in signalList: continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %15s'%'process')
        for i in xrange(0,len(signalList)) : datacard.write('%15s'%str(i+1-len(signalList)))
        i=0
        for proc in exp: 
            if proc in signalList: continue
            i=i+1
            datacard.write('%15s'%str(i))
        datacard.write('\n')
        datacard.write('\t\t\t %15s'%'rate')
        for sig in signalList: datacard.write('%15s'%('%3.2f'%exp[sig].Integral()))
        for proc in exp: 
            if proc in signalList: continue
            datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
        datacard.write('\n')
        datacard.write('-'*50+'\n')

        nomShapes=exp.copy()
        nomShapes['data_obs']=obs
        saveToShapesFile(outFile,nomShapes,'nom')
    
        #experimental systematics
        _,expVarShapes=getDistsFrom(directory=fIn.Get('%sshapes_%s_exp'%(opt.dist,cat)))
        nExpSysts=expVarShapes[signalList[0]].GetNbinsY()/2
        for isyst in xrange(1,nExpSysts):
            
            #test which variations are significant
            bwList={}
            ybin=2*(isyst-1)+1
            systVar=''
            upShapes,downShapes={},{}
            for proc in exp:
                if proc in expVarShapes:
                    
                    systVarDown = expVarShapes[proc].GetYaxis().GetBinLabel(ybin)
                    systVarUp   = expVarShapes[proc].GetYaxis().GetBinLabel(ybin+1)
                    systVar     = systVarUp[:-2]

                    downShapeH  = expVarShapes[proc].ProjectionX('%s%dDown'%(proc,isyst), ybin,   ybin)
                    upShapeH    = expVarShapes[proc].ProjectionX('%s%dUp'%(proc,isyst),   ybin+1, ybin+1)

                    bwList[proc]     = acceptVariationForDataCard(nomH=exp[proc], upH=upShapeH, downH=downShapeH)
                    if not bwList[proc]: continue

                    downShapes[proc] = downShapeH
                    upShapes[proc]   = upShapeH
                    
            #test if at least one process has been white listed
            if len(upShapes)+len(downShapes)==0:
                print 'Skipping',systVar,'for %s'%cat
                continue

            #export to shapes file
            saveToShapesFile(outFile,downShapes,systVar+'Down')
            saveToShapesFile(outFile,upShapes,systVar+'Up')

            #write to datacard
            datacard.write('%32s shapeN2'%systVar)        
            for sig in signalList: 
                if sig in bwList and bwList[sig]:
                    datacard.write('%15s'%'1') 
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if proc in bwList and bwList[proc] :
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        #rate systematics
        rateSysts=[
            ('lumi',           1.046,    'lnN',    []                ,['Multijetsdata']),
            ]
        try:
            jetCat=cat[:-2] if cat.endswith('t') else cat
            rateSysts.append( ('MultiJetsNorm%s'%jetCat,  qcdNormUnc[jetCat],  'lnU',    ['Multijetsdata']     ,[]) )
        except:
            pass

        for syst,val,pdf,whiteList,blackList in rateSysts:

            datacard.write('%32s %8s'%(syst,pdf))
            entryTxt=''
            try:
                entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
            except:
                entryTxt='%15s'%('%3.3f'%val)

            for sig in signalList: 
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        #generator level systematics 
        if systfIn is None: continue
        sampleSysts=[
            ('Mtop',            {'tbart'         : ['tbartm=169.5','tbartm=175.5'],  'tW':['tWm=169.5','tWm=175.5'] },                True , True),
            ('ttPartonShower',  {'tbart'         : ['tbartscaledown','tbartscaleup']},                                                True , True),
            ('tWscale',         {'tW'            : ['tWscaledown','tWscaleup']},                                                      True , True),            
            ('NLOgenerator',    {'tbart'         : ['tbartaMCNLO']},                                                                True , True),
            #('Hadronizer',      {'tbart'         : ['tbartaMCatNLO']},                                                                True , True),
            ('wFactScale',           { 'Wl': ['mur1muf0.5','mur1muf2'],   'Wc': ['mur1muf0.5','mur1muf2'],  'Wb': ['mur1muf0.5','mur1muf2'] },   False, False),
            ('wRenScale',            { 'Wl': ['mur0.5muf1','mur2muf1'],   'Wc': ['mur0.5muf1','mur2muf1'],  'Wb': ['mur0.5muf1','mur2muf1'] },   False, False),
            ('wCombScale',           { 'Wl': ['mur0.5muf0.5','mur2muf2'], 'Wc': ['mur0.5muf0.5','mur2muf2'],'Wb': ['mur0.5muf0.5','mur2muf2'] }, False, False),
            ('ttFactScale',          { 'tbart': ['muR1muF0.5','muR1muF2'] },                                                                     False , False),
            ('ttRenScale',           { 'tbart': ['muR0.5muF1','muR2muF1'] },                                                                     False , False),
            ('ttCombScale',          { 'tbart': ['muR0.5muF0.5','muR2muF2'] },                                                                    False , False),
            ]
        _,genVarShapes = getDistsFrom(directory=fIn.Get('%sshapes_%s_gen'%(opt.dist,cat)))
        _,altExp       = getDistsFrom(directory=systfIn.Get('%s_%s'%(opt.dist,cat)))
        for systVar, procsToApply, normalize, useAltShape in sampleSysts:

            #prepare shapes and check if variation is significant
            downShapes, upShapes = {}, {}
            
            for iproc in procsToApply:

                nomH=exp[iproc]

                #check which shape to use
                if useAltShape:

                    #get directly from another file
                    downH  = altExp[ procsToApply[iproc][0] ]
                    if len( procsToApply[iproc] ) > 1 :
                        upH    = altExp[ procsToApply[iproc][1] ]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
                else:

                    #project from 2D histo (re-weighted from nominal sample)
                    ybinUp, ybinDown = -1, -1
                    for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
                        label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
                        if procsToApply[iproc][0] in label : ybinDown=ybin
                        if procsToApply[iproc][1] in label : ybinUp=ybin

                    downH = genVarShapes[ iproc ].ProjectionX('%s%sDown'%(iproc,systVar), ybinDown, ybinDown)
                    upH   = genVarShapes[ iproc ].ProjectionX('%s%sUp'%(iproc,systVar),   ybinUp,   ybinUp)

                #normalize (shape only variation is considered)
                if normalize : downH.Scale( nomH.Integral()/downH.Integral() ) 
                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                #check if variation is meaningful
                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                if not accept : continue
                
                #save
                downShapes[iproc]=downH
                upShapes[iproc]=upH

            #check if something has been accepted
            if len(upShapes)==0 : continue

            #export to shapes file
            saveToShapesFile(outFile,downShapes,systVar+'Down')
            saveToShapesFile(outFile,upShapes,systVar+'Up')

            #write to datacard
            datacard.write('%32s shapeN2'%systVar)
            for sig in signalList: 
                if sig in procsToApply:
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if proc in procsToApply:
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')


        #all done
        datacard.close()




"""                                                                                                                                                                                                               
for execution from another script                                                                                                                                                                           
"""
if __name__ == "__main__":
    sys.exit(main())
