import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *

"""
Wrapper to be used when run in parallel
"""
def RunMethodPacked(args):
    inF,outF,channel,charge,wgt,isTT,flav,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr=args
    try:
        ROOT.ReadTree(str(inF),str(outF),channel,charge,wgt,isTT,flav,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    ROOT.AutoLibraryLoader.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    ROOT.gROOT.LoadMacro('src/ReadTree.cc+')
    from ROOT import ReadTree


    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='outDir',      help='output directory',                           default='analysis', type='string')
    parser.add_option('-j', '--json',        dest='json',        help='json file to process',                       default=None,       type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process',             default=None,       type='string')
    parser.add_option(      '--resetCache',  dest='resetCache',  help='reset normalization cache',                  default=False,       action='store_true')
    parser.add_option(      '--ch',          dest='channel',     help='channel',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge',                                     default=0,          type=int)
    parser.add_option(      '--flav',        dest='flav',        help='flavour splitting (for single files)',       default=0,          type=int)
    parser.add_option(      '--isTT',        dest='isTT',        help='ttbar sample (for single files)',            default=False,      action='store_true')
    parser.add_option(      '--norm',        dest='norm',        help='normalization scale (for single files)',     default=1.0,        type=float)
    parser.add_option(      '--genWgtMode',  dest='genWgtMode',  help='gen level wgts 0=none, 1=gen weight (for single files)',     default=0,        type=int)
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
    (opt, args) = parser.parse_args()

    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #prepare output
    os.system('mkdir -p %s'%opt.outDir)
    
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        if '/store/' in opt.input: opt.input='root://eoscms.cern.ch//eos/cms/'+opt.input
        outF=os.path.join(opt.outDir,os.path.basename(opt.input))
        task_list.append( (opt.input,outF,opt.channel,opt.charge,opt.norm,opt.isTT,opt.flav,opt.genWgtMode) )
    else:
        #read list of samples
        jsonFile = open(opt.json,'r')
        samplesList=json.load(jsonFile,encoding='utf-8').items()
        jsonFile.close()

        #read normalization
        xsecWgts, integLumi, puWgts = {}, {}, {}
        cache='%s/.xsecweights.pck'%opt.outDir
        if opt.resetCache : 
            print 'Removing normalization cache'
            os.system('rm %s'%cache)
        try:
            cachefile = open(cache, 'r')
            xsecWgts  = pickle.load(cachefile)
            integLumi = pickle.load(cachefile)
            puWgts    = pickle.load(cachefile)
            cachefile.close()        
            print 'Normalization read from cache (%s)' % cache
            for tag,sample in samplesList:
                if not tag in xsecWgts:
                    raise KeyError
        except:
            print 'Computing original number of events and storing in cache, this may take a while if it\'s the first time'
            xsecWgts,integLumi, puWgts = produceNormalizationCache(samplesList=samplesList,inDir=opt.input,cache=cache)
            
        #create the analysis jobs
        for tag,sample in samplesList:

            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if itag in tag:
                        processThisTag=True
                if not processThisTag : continue
            processedTags.append(tag)

            #get configuration
            isTT=sample[5]
            doFlavourSplitting=sample[6]
            genWgtMode=sample[7]
            wgt = xsecWgts[tag]
            puWgtGr,puUpWgtGr,puDownWgtGr=None,None,None
            try:
                puWgtGr,puUpWgtGr,puDownWgtGr=puWgts[tag]
            except:
                pass
            input_list=getEOSlslist(directory=opt.input+'/'+tag)
            for ifctr in xrange(0,len(input_list)):
                inF=input_list[ifctr]
                outF=os.path.join(opt.outDir,'%s_%d.root' % (tag,ifctr) )
                if doFlavourSplitting:
                    for flav in [0,1,4,5]:
                        task_list.append( (inF,outF,opt.channel,opt.charge,wgt,isTT,flav,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr) )                
                else:
                    task_list.append( (inF,outF,opt.channel,opt.charge,wgt,isTT,0,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr) )                

    #run the analysis jobs
    if opt.njobs == 0:
        for inF,outF,channel,charge,wgt,isTT,flav,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr in task_list:      
            ROOT.ReadTree(str(inF),str(outF),channel,charge,wgt,isTT,flav,genWgtMode,puWgtGr,puUpWgtGr,puDownWgtGr)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(RunMethodPacked, task_list)

    #merge final outputs
    for tag,sample in samplesList:
        if not tag in processedTags: continue
        doFlavourSplitting=sample[6]
        groupsToHadd=[tag]
        if doFlavourSplitting : 
            for flav in [1,4,5]:
                grouptsToHadd.append('%d_%s'%(flav,tag))
        for itag in groupsToHadd:
            os.system('hadd -f %s/%s.root %s/%s_*.root' % (opt.outDir,itag,opt.outDir,itag) )
            os.system('rm %s/%s_*.root' % (opt.outDir,itag) )
    print 'Analysis results are available in %s' % opt.outDir


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
