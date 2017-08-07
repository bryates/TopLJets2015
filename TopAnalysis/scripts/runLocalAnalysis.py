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


    method,inF,outF,channel,charge,flav,runSysts,era,runPeriod,tag,debug=args
    print 'Running ',method,' on ',inF
    print 'Output file',outF
    print 'Selection ch=',channel,' charge=',charge,' flavSplit=',flav,' systs=',runSysts
    print 'Normalization applied from tag=',tag
    print 'Corrections will be retrieved for era=',era
    print 'Run period',runPeriod
    if debug : print 'Verbose mode'

    try:
        cmd='analysisWrapper --era %s --runPeriod %s --normTag %s --in %s --out %s --method %s --charge %d --channel %d --flav %d' %(era,
                                                                                                                      runPeriod,
                                                                                                                      tag,
                                                                                                                      inF,
                                                                                                                      outF,
                                                                                                                      method,
                                                                                                                      charge,
                                                                                                                      channel,
                                                                                                                      flav)
        if runSysts : cmd += ' --runSysts'
        if debug : cmd += ' --verbose'
        print cmd
        os.system(cmd)
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
    parser.add_option('-m', '--method',      dest='method',      help='method to run [%default]',                   default='TOP-16-006::RunTop16006',  type='string')
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file [%default]',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process  [%default]',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--flav',        dest='flav',        help='split according to heavy flavour content  [%default]',   default=0,          type=int)
    parser.add_option(      '--ch',          dest='channel',     help='channel  [%default]',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge  [%default]',                                     default=0,          type=int)
    parser.add_option(      '--era',         dest='era',         help='era to use for corrections/uncertainties  [%default]',   default='era2016',  type='string')
    parser.add_option(      '--runPeriod',   dest='runPeriod',   help='peirod to use for corrections/uncertainties  [%default]',default='BCDEF',    type='string')
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag  [%default]',                    default=None,       type='string')
    parser.add_option(      '--type',        dest='type',        help='only run Data or MC  [%default]',                        default=None,       type='string')
    parser.add_option(      '--dataOnly',    dest='dataOnly',    help='only run Data  [%default]',                              default=False,      action='store_true')
    parser.add_option(      '--MCOnly',      dest='MCOnly',      help='only run MC  [%default]',                                default=False,      action='store_true')
    parser.add_option('-q', '--queue',       dest='queue',       help='submit to this queue  [%default]',                       default='local',    type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                  default=0,          type='int')
    parser.add_option('-v', '--verbose',     dest='debug',       help='pint debug messages [%default]',                         default=False,      action='store_true')
    (opt, args) = parser.parse_args()

    #parse selection list
    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #prepare output if a directory
    if not '.root' in opt.output:
        os.system('mkdir -p %s'%opt.output)
    
    #correct location of corrections to be used using cmsswBase, if needed
    cmsswBase=os.environ['CMSSW_BASE']
    if not cmsswBase in opt.era : opt.era=cmsswBase+'/src/TopLJets2015/TopAnalysis/data/'+opt.era
    era=opt.era

    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : inF='root://eoscms//eos/cms'+opt.input        
        print inF
        outF=opt.output
        task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flav,opt.runSysts,opt.era,opt.runPeriod,opt.tag,opt.debug) )
    else:

        inputTags=getEOSlslist(directory=opt.input,prepend='')
        for baseDir in inputTags:

            tag=os.path.basename(baseDir)
            if tag=='backup' : continue

            #filter tags
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if itag[0]=='^': 
                        itag=itag[1:]
                        if itag in tag : 
                            processThisTag=False
                            break
                        else : 
                            processThisTag=True
                    elif itag in tag:
                        processThisTag=True
                        break
                if not processThisTag : continue
            splitRun = lambda x: ["2016" + x[i] for i in range(0, len(x), 1)]
            split_run = splitRun( opt.runPeriod )
            if 'Data' in tag:
                if not any(run in tag for run in split_run): continue
                if opt.MCOnly: continue
            if 'MC' in tag and opt.dataOnly: continue

            ############### Submit to condor ###############
            input_list=getEOSlslist(directory='%s/%s' % (opt.input,tag) )
            ext_len=",".join(input_list).count("_ext")
            target = '%s/src/TopLJets2015/TopAnalysis/%s' % (cmsswBase,tag)
            condorFile = open(target,'w')
            condorFile.write('universe              = vanilla\n')
            condorFile.write('executable            = condor/cond_submit.sh\n')
            condorFile.write('arguments             = $(ClusterID) $(ProcId) %s %s %s %s\n' % (opt.input,opt.output,tag,opt.runPeriod))
            condorFile.write('output                = condor/log/%s_$(ProcId).out\n' % tag)
            condorFile.write('error                 = condor/log/%s_$(ProcId).err\n' % tag)
            condorFile.write('log                   = condor/log/%s.log\n' % tag)
            condorFile.write('+JobFlavour           = "workday"\n')
            condorFile.write('Should_Transfer_Files = NO\n')
            condorFile.write('queue %d' % (len(input_list)-ext_len))
            condorFile.close()
            os.system('condor_submit %s -batch-name %s' % (target,tag))
            os.system('rm %s' % (tag))
            ############### Special case for ext samples ###############
            if "_ext" in tag:
              target = '%s/src/TopLJets2015/TopAnalysis/%s' % (cmsswBase,tag)
              condorFile = open(target,'w')
              condorFile.write('universe              = vanilla\n')
              condorFile.write('executable            = condor/cond_ext.sh\n')
              condorFile.write('arguments             = $(ClusterID) $(ProcId) %s %s %s %s\n' % (opt.input,opt.output,tag,opt.runPeriod))
              condorFile.write('output                = condor/log/%s_$(ProcId).out\n' % tag+"_ext")
              condorFile.write('error                 = condor/log/%s_$(ProcId).err\n' % tag+"_ext")
              condorFile.write('log                   = condor/log/%s.log\n' % tag+"_ext")
              condorFile.write('+JobFlavour           = "workday"\n')
              condorFile.write('Should_Transfer_Files = NO\n')
              condorFile.write('queue %d' % ext_len)
              condorFile.close()
              os.system('condor_submit %s -batch-name %s' % (target,tag))
              os.system('rm %s' % (tag))
            for ifile in xrange(0,len(input_list)):
                inF=input_list[ifile]
                outF=os.path.join(opt.output,'%s_%d.root' %(tag,ifile))
                task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flav,opt.runSysts,opt.era,opt.runPeriod,tag,opt.debug) )

    #run the analysis jobs
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for args in task_list: RunMethodPacked(args)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        queue = opt.queue
        if queue is "8nh": queue = "longlunch"
        print 'launching %d tasks to submit to the %s queue'%(len(task_list),queue)
        for method,inF,outF,channel,charge,flav,runSysts,era,runPeriod,tag,debug in task_list:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py -i %s -o %s --charge %d --ch %d --era %s --runPeriod %s --tag %s --flav %d --method %s' % (cmsswBase,inF,outF,charge,channel,era,runPeriod,tag,flav,method)
            if debug : localrun += ' --verbose %s' % (debug)
            if runSysts : localRun += ' --runSysts'            
            ############### Now using condor instead of LSF (bsub) ###############
            #cmd='bsub -q %s -J %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,outF,cmsswBase,localRun)
            #os.system(cmd)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
