import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr

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
in case of multiple signals, remove others
"""
def filterShapeList(exp,signalList,rawSignalList):
    newExp={}
    for key in exp:

        matchFound=False
        for rawKey in rawSignalList:
            if rawKey in key:
                matchFound=True
        if matchFound and not key in signalList : continue

        newExp[key]=exp[key]
    return newExp


"""
get distributions from file
"""
def getDistsFrom(directory,addList=None,addTitle='',keyFilter=''):
    obs=None
    exp={}
    dirName=directory.GetName()
    addedHistProduced=False
    for key in directory.GetListOfKeys():
        if len(keyFilter)>0 and key.GetName()!='%s_%s'%(dirName,keyFilter) : continue
        obj=directory.Get(key.GetName())
        if not obj.InheritsFrom('TH1') : continue
        inAddList = False;
        if obj.GetName()==dirName :
            obs=obj.Clone('data_obs')
            obs.SetDirectory(0)
        else :
            newName=obj.GetName().split(dirName+'_')[-1]
            for token in ['+','-','*',' ','#','{','(',')','}','@']:
                newName=newName.replace(token,'')
            if addList and newName in addList :
                newName='added' if addTitle == '' else addTitle
                if addedHistProduced :
                    exp[newName].Add(obj)
                    addedHistProduced = True;
                    inAddList=True
                else :
                    exp[newName] = obj.Clone(newName)
            else :
                exp[newName]=obj.Clone(newName)
                exp[newName].SetDirectory(0)

            #newName = ('added' if addTitle == '' else addTitle) if inAddList and addedHistProduced else newName
            if obj.InheritsFrom('TH2'):
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
merge many directories with the same directory base-name
  '.' --> 'p'
"""
def getMergedDists(fIn,basedir='',addList=None,addTitle='',
                   widList=None,nomWid='',sigList=None,keyFilter=''):
    obs=None
    exp={}
    for wid in widList :
        dirName='%s%s'%(basedir,wid)
        print dirName
        directory = fIn.Get(dirName)
        if wid == nomWid :
            obs,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                newName=key
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                print " - " +newName
                exp[newName]=mergeExp[key].Clone(newName)
        else :
            _,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                    print " - " +newName
                    exp[newName]=mergeExp[key].Clone(newName)
                else : continue
    return obs,exp


"""
save distributions to file
"""
def saveToShapesFile(outFile,shapeColl,directory=''):
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    if len(directory)==0:
        fOut.cd()
    else:
        fOut.mkdir(directory)
        outDir=fOut.Get(directory)
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
    parser.add_option('-i', '--input',     dest='input',     help='input plotter',                                    default=None,            type='string')
    parser.add_option(      '--systInput', dest='systInput', help='input plotter for systs from alternative samples', default=None,            type='string')
    parser.add_option('-d', '--dist',      dest='dist',      help='distribution',                                     default='mlb',           type='string')
    parser.add_option('-s', '--signal',    dest='signal',    help='signal (csv)',                                     default='tbart,tW',      type='string')
    parser.add_option('-w', '--wids',      dest='widList',   help='signal widths',                                    default='0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0',type='string')
    parser.add_option('-c', '--cat',       dest='cat',       help='categories (csv)',                                 default='1b,2b',         type='string')
    parser.add_option('-o', '--output',    dest='output',    help='output directory',                                 default='datacards',     type='string')
    parser.add_option(      '--addSigs',   dest='addSigs',   help='signal processes to add',                          default=False,           action='store_true')
    parser.add_option(      '--lfs',       dest='lfsInput',  help='lepton final states to consider',                  default='EE,EM,MM',  type='string')
    parser.add_option(      '--lbCat',     dest='lbCat',     help='lepton final states to consider',                  default='highpt,lowpt',  type='string')
    (opt, args) = parser.parse_args()

    # parse the channels, lb categories to consider
    lfsList=opt.lfsInput.split(',')
    lbCatList=opt.lbCat.split(',')

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')
    signalList = { 'top' } if opt.addSigs else rawSignalList[:]

    # what are our widths?
    rawWidList = opt.widList.split(',')
    widList = rawWidList[:]
    for i in xrange(0,len(rawWidList)) :
        widList[i] = ('%sw'%(rawWidList[i]))
    nomWid="1.0w"
    modNomWid="1p0w"

    modWidList=widList[:]
    for i in xrange(0,len(widList)) :
        modWidList[i]=widList[i].replace('.','p')

    # what are our categories?
    catList=opt.cat.split(',')

    # prepare systematics
    tWRateSystList=['tW']
    ttRateSystList=['tbart']
    tWIsSig=('tW'    in rawSignalList)
    ttIsSig=('tbart' in rawSignalList)
    for wid in modWidList :
        if tWIsSig : tWRateSystList += ["tW%s"%wid]
        if ttIsSig : ttRateSystList += ["tbart%s"%wid]

    rateSysts=[
          ('lumi_13TeV',       1.027,    'lnN',    []                   ,['Multijetsdata']),
          #('DYnorm_th',        1.038,    'lnN',    ['DYl','DYc','DYb']  ,[]),
          #('Wnorm_th',         1.037,    'lnN',    ['Wl' ,'Wc','Wb']    ,[]),
          ('DYnorm_th',        1.038,    'lnN',    ['DY']  ,[]),
          ('Wnorm_th',         1.037,    'lnN',    ['W']   ,[]),
          ('tWnorm_th',        1.054,    'lnN',    tWRateSystList,[]),
          ('tnorm_th',         1.044,    'lnN',    ['tch']              ,[]),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,[]),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,[]),
    ]

    Mtop=['tbartm=169.5','tbartm=175.5','tWm=169.5','tWm=175.5']
    ttParton_tt=['tbartscaledown','tbartscaleup']
    ttParton_tW=['tWscaledown','tWscaleup']
    tWinterf=['tWDS']

    MtopMap={}
    ttPartonMap={}
    tWinterfMap={}
    for wid in modWidList :
        MtopMap[    'tbart%s'%wid]=['tbartm=169.5%s'  %(wid),'tbartm=175.5%s'%(wid)]
        ttPartonMap['tbart%s'%wid]=['tbartscaledown%s'%(wid),'tbartscaleup%s'%(wid)]
        MtopMap[    'tW%s'   %wid]=['tWm=169.5%s'     %(wid),'tWm=175.5%s'   %(wid)]
        ttPartonMap['tW%s'   %wid]=['tWscaledown%s'   %(wid),'tWscaleup%s'   %(wid)]

        tWinterfMap['tW%s'   %wid]=['tWDS%s'   %(wid)]

    sampleSysts=[
          #ttbar modelling
          ('Mtop'          , MtopMap    , True ,  True, False),
          ('ttPartonShower', ttPartonMap, False , True, False),
          #tWinterference
          ('tWttinterf'    , tWinterfMap, False , True, True)
    ]

    genSysts=[
        ('jes',  False,False,False),
        ('les',  False,False,False),
        ('jer',  False,False,False),
        ('btag', False,False,False),
        ('pu',   False,False,False)
    ]

    #prepare output directory
    os.system('mkdir -p %s'%opt.output)

    anCat=''
    for subDir in opt.input.split('/'):
        if 'analysis_' not in subDir: continue
        anCat=subDir.replace('analysis_','')

    #get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input)
    systfIn=None
    if opt.systInput:
        systfIn=ROOT.TFile.Open(opt.systInput)

    #prepare output ROOT file
    outFile='%s/shapes.root'%(opt.output)
    fOut=ROOT.TFile.Open(outFile,'RECREATE')
    fOut.Close()

    #loop over lepton final states
    for (lbCat,lfs,cat) in [(a,b,c) for a in lbCatList for b in lfsList for c in catList] :
        obs=ROOT.TH1F('','',100,0,100)
        exp={}
        #nominal expectations
        if opt.addSigs :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,opt.dist)),rawSignalList,'top',widList,nomWid,signalList)
        else :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,opt.dist)),None,'',widList,nomWid,signalList)
        #exp=filterShapeList(exp,signalList,rawSignalList) ???

        nomShapes=exp.copy()
        nomShapes['data_obs']=obs
        saveToShapesFile(outFile,nomShapes,('%s%s%s_%s'%(lbCat,lfs,cat,opt.dist)))

        #loop over categories, widths
        for wid in modWidList:
            if modNomWid in wid : continue

            print 'Initiating %s datacard for %s%s%s_%s'%(opt.dist,lbCat,lfs,cat,wid)


            #start the datacard
            datacard=open('%s/datacard__%s_%s%s%s.dat'%(opt.output,wid,lbCat,lfs,cat),'w')
            datacard.write('#\n')
            datacard.write('# Generated by %s with git hash %s for analysis category %s%s%s_%s\n' % (getpass.getuser(),
                commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                lbCat, lfs, cat, wid) )

            datacard.write('#\n')
            datacard.write('imax *\n')
            datacard.write('jmax *\n')
            datacard.write('kmax *\n')
            datacard.write('-'*50+'\n')
            datacard.write('shapes *        * shapes.root %s%s%s_mlb/$PROCESS $SYSTEMATIC/$PROCESS\n'%(lbCat,lfs,cat))
            datacard.write('-'*50+'\n')
            datacard.write('bin 1\n')
            datacard.write('observation %3.1f\n' % obs.Integral())
            datacard.write('-'*50+'\n')

            #expectations
            datacard.write('\t\t\t %32s'%'bin')
            for i in xrange(0,len(exp)+(2 if opt.addSigs else 0)-(len(modWidList)-1)*len(signalList)+len(signalList)):
                datacard.write('%15s'%'1')
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,modNomWid)))
            for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,wid)))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc : isSig=True
                if isSig : continue
                datacard.write('%15s'%proc)
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for i in xrange(0,2*len(signalList)) : datacard.write('%15s'%str(i+1-2*len(signalList)))
            i=0
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc : isSig=True
                if isSig : continue
                i=i+1
                datacard.write('%15s'%str(i))
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'rate')
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].Integral())))
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,wid)].Integral())))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc : isSig=True
                if isSig : continue
                datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
            datacard.write('\n')
            datacard.write('-'*50+'\n')


            #rate systematics
            try:
                jetCat=cat[:-2] if cat.endswith('t') else cat
                rateSysts.append( ('MultiJetsNorm%s%s'%(jetCat,anCat),  1+qcdNorm[jetCat][1],                       'lnN',    ['Multijetsdata']    ,[]) )
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
                    newSig='%s%s'%(sig,modNomWid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    newSig='%s%s'%(sig,wid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc : isSig=True
                    if isSig : continue
                    if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

            #generator level systematics
            if systfIn is None :
                datacard.close()
                continue

            systSignalList=Mtop+ttParton_tt+ttParton_tW+tWinterf

            # sample systematics
            _,genVarShapes = getMergedDists(fIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,opt.dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)
            _,altExp       = getMergedDists(systfIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,opt.dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)

            for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:

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

                    # use do down/up x nom to generate the variation, then mirror it
                    if projectRelToNom:
                        ratioH=downH.Clone()
                        ratioH.Divide(upH)
                        for xbin in xrange(1,nomH.GetNbinsX()+1):
                            nomVal=nomH.GetBinContent(xbin)
                            varVal = ratioH.GetBinContent(xbin) * nomVal
                            upH.SetBinContent(xbin, varVal)
                            varVal = varVal- nomVal
                            downH.SetBinContent(xbin, nomVal-varVal)

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
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc : isSig=True
                    if isSig : continue
                    if proc in procsToApply and proc in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')


            #
            for systVar, normalize, useAltShape, projectRelToNom in genSysts:
                _,genVarShapesUp = getMergedDists((systfIn if useAltShape else fIn),
                        ('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,opt.dist)),
                        (rawSignalList if opt.addSigs else None),
                        ('top' if opt.addSigs else ''),
                        widList,nomWid,signalList)
                _,genVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                        ('%sdn%s%s%s_%s_'%(systVar,lbCat,lfs,cat,opt.dist)),
                        (rawSignalList if opt.addSigs else None),
                        ('top' if opt.addSigs else ''),
                        widList,nomWid,signalList)

                #prepare shapes and check if variation is significant
                downShapes, upShapes = {}, {}

                for iproc in exp:
                    nomH=exp[iproc]

                    if not iproc in genVarShapesDn :
                        continue
                    #get directly from another file
                    downH  = genVarShapesDn[iproc]
                    if iproc in genVarShapesUp:
                        upH = genVarShapesUp[iproc]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)

                    # use do down/up x nom to generate the variation, then mirror it
                    if projectRelToNom:
                        ratioH=downH.Clone()
                        ratioH.Divide(upH)
                        for xbin in xrange(1,nomH.GetNbinsX()+1):
                            nomVal=nomH.GetBinContent(xbin)
                            varVal = ratioH.GetBinContent(xbin) * nomVal
                            upH.SetBinContent(xbin, varVal)
                            varVal = varVal- nomVal
                            downH.SetBinContent(xbin, nomVal-varVal)

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
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in exp and ("%s%s"%(sig,wid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc : isSig=True
                    if isSig : continue
                    if proc in exp and proc in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

        #all done
        datacard.close()

    ######
    # BEGIN MORPHING
    ######

    if systfIn is None : return

    print "\n Creating morphed dists\n"

    widths=map(float,rawWidList)
    masses=[(169.5,opt.systInput,'t#bar{t} m=169.5'),
            (172.5,opt.input    ,'t#bar{t}'),
            (175.5,opt.systInput,'t#bar{t} m=175.5'),
            (169.5,opt.systInput,'tW m=169.5'),
            (172.5,opt.input    ,'tW'),
            (175.5,opt.systInput,'tW m=175.5')]

    # make workspace for each signal process
    for sig in signalList :
        modSig=replaceBadCharacters(sig)

        ws=ROOT.RooWorkspace('sigws_%s'%modSig)
        ws.factory('x[0,300]')

        widthDim=ROOT.RooBinning(len(widths)-1,0,len(widths)-1)
        massDim=ROOT.RooBinning(len(masses)-1,0,len(masses)-1)
        refGrid=ROOT.RooMomentMorphND.Grid(widthDim,massDim)

        for lbCat,ch,cat in [(lbCat,ch,cat) for lbCat in ['highpt'] for ch in ['EM'] for cat in ['2b']]:
            for imass in xrange(0,len(masses)):
                mass,url,proc=masses[imass]

                if modSig not in replaceBadCharacters(proc) : continue

                tfIn=ROOT.TFile.Open(url)
                for iwid in xrange(0,len(widths)):

                    #get histogram and convert to a PDF
                    dirname='%s%s%s_%s_%3.1fw'%(lbCat,ch,cat,opt.dist,widths[iwid])
                    h=tfIn.Get(dirname+"/"+dirname+'_'+proc)
                    name='%s_m%d'%(dirname,int(10*mass))
                    data=ROOT.RooDataHist(name,name,ROOT.RooArgList(ws.var("x")),h)
                    pdf=ROOT.RooHistPdf(name+"_pdf",name+"_pdf",ROOT.RooArgSet(ws.var("x")),data)
                    getattr(ws,'import')(pdf,ROOT.RooCmdArg())

                    #add pdf to the grid
                    print 'Adding',pdf.GetName(),'@ (',iwid,imass,')'
                    refGrid.addPdf(ws.pdf(pdf.GetName()),iwid,imass)
                tfIn.Close()

        # produce morphed pdf
        ws.factory('alpha[0,3]')
        ws.factory('beta[0,3]')
        pdf=ROOT.RooMomentMorphND('widmorphpdf','widmorphpdf',
                                  ROOT.RooArgList( ws.var('alpha'), ws.var('beta') ),
                                  ROOT.RooArgList( ws.var('x') ),
                                  refGrid,
                                  ROOT.RooMomentMorphND.Linear)
        pdf.useHorizontalMorphing(False)
        getattr(ws,'import')(pdf,ROOT.RooCmdArg())

        ws.var('alpha').setVal(1.0)
        ws.var('beta').setVal(1.0)

        # save workspace to shapes
        outFile='%s/shapes.root'%(opt.output)
        fOut=ROOT.TFile.Open(outFile,'UPDATE')
        fOut.cd()
        ws.Write()
        fOut.Close()

    print "\n Morphed dists created, workspace saved. Producing datacards: \n"

    #####
    # START WIDTH/MASS SCAN DATACARD
    #####
    fIn=ROOT.TFile.Open(opt.input)
    for lbCat,ch,cat in [(lbCat,ch,cat) for lbCat in lbCatList for ch in lfsList for cat in catList] :
        # get signal/backgrounds
        obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,opt.dist)),None,'',widList,nomWid,signalList)

        # prepare datacard
        dirname='%s%s%s_mlb_%3.1fw'%(lbCat,ch,cat,widths[iwid])
        datacard=open('%s/datacard_widfit__%s%s%s.dat'%(opt.output,lbCat,ch,cat),'w')

        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis category %s%s\n' % (getpass.getuser(),
            commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
            cat, ch) )

        # header
        datacard.write('#\n')
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')

        # shapes locations
        datacard.write('-'*50+'\n')
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc : isSig=True
            if isSig : continue
            datacard.write('shapes %10s * shapes.root %s%s%s_mlb/$PROCESS $SYSTEMATIC/$PROCESS\n'%(proc,lbCat,ch,cat))

        datacard.write('shapes data_obs   * shapes.root %s%s%s_mlb/$PROCESS $SYSTEMATIC/$PROCESS\n'%(lbCat,ch,cat))

        for sig in signalList:
            datacard.write('shapes %10s * shapes.root sigws_$PROCESS:widmorphpdf\n'%(sig))
        datacard.write('-'*50+'\n')

        # bin names
        datacard.write('bin 1\n')
        datacard.write('observation %3.1f\n' % obs.Integral())

        # observations
        datacard.write('-'*50+'\n')
        datacard.write('\t\t\t %32s'%'bin')
        for i in xrange(0,len(exp)-(len(widths)-1)*len(signalList)): datacard.write('%15s'%'1')
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'process')
        for sig in signalList:
            datacard.write('%15s'%(sig))
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc : isSig=True
            if isSig : continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'process')
        for i in xrange(0,len(signalList)) : datacard.write('%15s'%str(-1*(len(signalList)-1-i)))
        i=0
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc : isSig=True
            if isSig : continue
            i=i+1
            datacard.write('%15s'%str(i))
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'rate')
        for sig in signalList :
            datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].Integral())))
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc : isSig=True
            if isSig : continue
            datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
        datacard.write('\n')
        datacard.write('-'*50+'\n')

        datacard.close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
