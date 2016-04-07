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
    inF,outF,channel,charge,wgtH,flav,runSysts=args
    try:
        ROOT.RunTop16006(str(inF),str(outF),channel,charge,flav,wgtH,runSysts)
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
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics',                            default=False,      action='store_true')
    parser.add_option(      '--cache',       dest='cache',       help='use this cache file',                        default='data/genweights.pck', type='string')
    parser.add_option(      '--flav',        dest='flav',        help='split according to heavy flavour content',   default=0,          type=int)
    parser.add_option(      '--ch',          dest='channel',     help='channel',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge',                                     default=0,          type=int)
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='submit to this queue',                       default='local',    type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',                                default=0,    type='int')
    (opt, args) = parser.parse_args()

    #compile macro
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    ROOT.gROOT.LoadMacro('src/TOP-16-006.cc+')
    from ROOT import RunTop16006

    #parse selection list
    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #prepare output if a directory
    if not '.root' in opt.output:
        os.system('mkdir -p %s'%opt.output)

    #read normalization
    cachefile = open(opt.cache, 'r')
    genWgts   = pickle.load(cachefile)
    cachefile.close()        
    print 'Normalization read from cache (%s)' % opt.cache
    
    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : inF='root://eoscms//eos/cms'+opt.input        
        outF=opt.output
        wgt=None
        if opt.tag :
            if opt.tag in genWgts:
                wgtH=genWgts[opt.tag]
        print inF,outF,opt.channel,opt.charge,wgtH,opt.flav,opt.runSysts
        task_list.append( (inF,outF,opt.channel,opt.charge,wgtH,opt.flav,opt.runSysts) )
    else:

        inputTags=getEOSlslist(directory=opt.input,prepend='')
        for baseDir in inputTags:

            tag=os.path.basename(baseDir)
            if tag=='backup' : continue

            #filter tags
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if itag in tag:
                        processThisTag=True
                if not processThisTag : continue

            wgtH=genWgts[tag] if opt.queue=='local' else tag
            input_list=getEOSlslist(directory='%s/%s' % (opt.input,tag) )
            for ifile in xrange(0,len(input_list)):
                inF=input_list[ifile]
                outF=os.path.join(opt.output,'%s_%d.root' %(tag,ifile))
                doFlavourSplitting=True if ('MC13TeV_WJets' in tag or 'MC13TeV_DY50toInf' in tag) else False
                if doFlavourSplitting:
                    for flav in [0,1,4,5]:
                        task_list.append( (inF,outF,opt.channel,opt.charge,wgtH,flav,opt.runSysts) )
                else:
                    task_list.append( (inF,outF,opt.channel,opt.charge,wgtH,0,opt.runSysts) )

    #run the analysis jobs
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for inF,outF,channel,charge,wgtH,flav,runSysts in task_list:
                ROOT.RunTop16006(str(inF),str(outF),channel,charge,flav,wgtH,runSysts)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        print 'launching %d tasks to submit to the %s queue'%(len(task_list),opt.queue)
        cmsswBase=os.environ['CMSSW_BASE']
        for inF,outF,channel,charge,tag,flav,runSysts in task_list:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py -i %s -o %s --charge %d --ch %d --tag %s --flav %d' % (cmsswBase,inF,outF,charge,channel,tag,flav)
            if runSysts : localRun += ' --runSysts'            
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            os.system(cmd)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
