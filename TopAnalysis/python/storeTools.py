import ROOT
import pickle

"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getEOSlslist(directory, mask='', prepend='root://eoscms//eos/cms'):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    eos_cmd = '/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select'
    data = Popen([eos_cmd, 'ls', '/eos/cms/'+directory],stdout=PIPE)
    out,err = data.communicate()

    full_list = []

    ## if input file was single root file:
    if directory.endswith('.root'):
        if len(out.split('\n')[0]) > 0:
            return [prepend + directory]

    ## instead of only the file name append the string to open the file in ROOT
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        full_list.append(prepend + directory + '/' + line)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        return stripped_list

    ## return 
    return full_list

"""
Loops over a list of samples and produces a cache file to normalize MC
"""
def produceNormalizationCache(samplesList,inDir,cache):

    xsecWgts,integLumi={},{}
    for tag,sample in samplesList: 

        if sample[1]==1 : 
            xsecWgts[tag]=1.0
            continue

        input_list=getEOSlslist(directory=inDir+'/'+tag)            
        xsec=sample[0]            
        norigEvents=0
        for f in input_list:
            fIn=ROOT.TFile.Open(f)
            norigEvents+=fIn.Get('analysis/counter').GetBinContent(1)
            fIn.Close()
        xsecWgts[tag]  = xsec/norigEvents  if norigEvents>0 else 0
        integLumi[tag] = norigEvents/xsec  if norigEvents>0 else 0
        print '... %s cross section=%f pb #orig events=%d lumi=%3.2f/fb' % (tag,xsec,norigEvents,integLumi[tag]/1000.)

        #dump to file    
        cachefile=open(cache,'w')
        pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(integLumi, cachefile, pickle.HIGHEST_PROTOCOL)
        cachefile.close()
        print 'Produced normalization cache (%s)'%cache
    return xsecWgts,integLumi
