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
    inF,outF,channel,charge,wgt,isTT,flav,genWgtMode,runSysts=args
    try:
        ROOT.ReadTree(str(inF),str(outF),channel,charge,wgt,isTT,flav,genWgtMode,runSysts)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='outDir',      help='output directory',                           default='analysis', type='string')
    parser.add_option('-j', '--json',        dest='json',        help='json file to process',                       default=None,       type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process',             default=None,       type='string')
    parser.add_option(      '--nocompile',   dest='nocompile',   help='don\'t compile macro',                       default=False,      action='store_true')
    parser.add_option(      '--resetCache',  dest='resetCache',  help='reset normalization cache',                  default=False,      action='store_true')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics',                            default=False,      action='store_true')
    parser.add_option(      '--cache',       dest='cache',       help='use this cache file',                        default='data/.xsecweights.pck', type='string')
    parser.add_option(      '--ch',          dest='channel',     help='channel',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge',                                     default=0,          type=int)
    parser.add_option(      '--flav',        dest='flav',        help='flavour splitting (for single files)',       default=0,          type=int)
    parser.add_option(      '--isTT',        dest='isTT',        help='ttbar sample (for single files)',            default=False,      action='store_true')
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='submit to this queue',                       default='local',    type='string')
    parser.add_option(      '--genWgtMode',  dest='genWgtMode',  help='gen level wgts 0=none, 1=gen weight (for single files)',   default=0,    type=int)
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',                                default=0,    type='int')
    (opt, args) = parser.parse_args()

    ROOT.AutoLibraryLoader.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    if opt.nocompile:
        ROOT.gROOT.LoadMacro('src/ReadTree.cc')
    else:
        ROOT.gROOT.LoadMacro('src/ReadTree.cc+')

    from ROOT import ReadTree



    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #prepare output
    os.system('mkdir -p %s'%opt.outDir)

    #read normalization
    xsecWgts, integLumi = {}, {}
    if opt.resetCache : 
        print 'Removing normalization cache'
        os.system('rm %s'%opt.cache)
    try:
        cachefile = open(opt.cache, 'r')
        xsecWgts  = pickle.load(cachefile)
        integLumi = pickle.load(cachefile)
        cachefile.close()        
        print 'Normalization read from cache (%s)' % opt.cache
    except:
        print 'Failed to read cache %s'%opt.cache


    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        if '/store/' in opt.input: opt.input='root://eoscms//eos/cms'+opt.input
        outF=os.path.join(opt.outDir,os.path.basename(opt.input))
        wgt=None
        if opt.tag :
            if opt.tag in xsecWgts:
                wgt=xsecWgts[opt.tag]
        task_list.append( (opt.input,outF,opt.channel,opt.charge,wgt,opt.isTT,opt.flav,opt.genWgtMode,opt.runSysts) )
    else:

        #read list of samples
        jsonFile = open(opt.json,'r')
        samplesList=json.load(jsonFile,encoding='utf-8').items()
        jsonFile.close()

        #check if all samples are in cache or force it's recreation
        try:

            if opt.resetCache : 
                raise KeyError

            for tag,sample in samplesList:
                if not tag in xsecWgts:
                    raise KeyError
        except:
            print 'Computing original number of events and storing in cache, this may take a while if it\'s the first time'
            print 'Starting with %d processes'%len(xsecWgts)
            xsecWgts, integLumi = produceNormalizationCache(samplesList=samplesList,
                                                            inDir=opt.input,
                                                            cache=opt.cache, 
                                                            xsecWgts=xsecWgts, 
                                                            integLumi=integLumi)
            print 'Current map has now %d processes'%len(xsecWgts)


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
            wgt = xsecWgts[tag] if opt.queue=='local' else tag
            input_list=getEOSlslist(directory=opt.input+'/'+tag)
            for ifctr in xrange(0,len(input_list)):
                inF=input_list[ifctr]
                outF=os.path.join(opt.outDir,'%s_%d.root' % (tag,ifctr) )
                if doFlavourSplitting:
                    for flav in [0,1,4,5]:
                        task_list.append( (inF,outF,opt.channel,opt.charge,wgt,isTT,flav,genWgtMode,opt.runSysts) )
                else:
                    task_list.append( (inF,outF,opt.channel,opt.charge,wgt,isTT,0,genWgtMode,opt.runSysts) )
                
    #run the analysis jobs
    if opt.queue=='local':
        if opt.njobs == 0:
            for inF,outF,channel,charge,wgt,isTT,flav,genWgtMode,runSysts in task_list:
                ROOT.ReadTree(str(inF),str(outF),channel,charge,wgt,isTT,flav,genWgtMode,runSysts)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        for inF,outF,channel,charge,wgt,isTT,flav,genWgtMode,runSysts in task_list:
            localOpt='--nocompile -i %s -o $s --charge %d --tag %s --isTT %d --flav %d --genWgtMode %d --runSysts %d' % (inF,outF,charge,wgt,isTT,flav,genWgtMode,runSysts)
            print localOpt

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
