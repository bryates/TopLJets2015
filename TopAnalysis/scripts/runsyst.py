import os

syst=["TRIGGER" ,"LEP" ,"PU" ,"PI" ,"JER"]
#syst=["TRIGGER" ,"LEP" ,"TRK" ,"PU" ,"PI" ,"JER"]
for i in xrange(0,len(syst)):
    cmd="python scripts/runLocalAnalysis.py -i /store/group/phys_top/byates/LJets2016/8db9ad6/ -o LJets2015/2016/etaPiK/up/%s --method TOP::RunTopKalman --era era2016 --runPeriod BCDEFGH --runSysts %d -q 8nh --MCOnly" % (syst[i],i+1)#2**i)
    print cmd
    os.system(cmd)
    cmd="python scripts/runLocalAnalysis.py -i /store/group/phys_top/byates/LJets2016/8db9ad6/ -o LJets2015/2016/etaPiK/down/%s --method TOP::RunTopKalman --era era2016 --runPeriod BCDEFGH --runSysts -%d -q 8nh --MCOnly" % (syst[i],i+1)#2**i)
    print cmd
    os.system(cmd)
