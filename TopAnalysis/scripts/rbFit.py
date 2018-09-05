import json
import math
from collections import OrderedDict

jsonFile = open('data/era2016/rbFit.json','r')
rbList=json.load(jsonFile,encoding='utf-8').items()#,object_pairs_hook=OrderedDict).items()
jsonFile.close()

nom=0
fsrdown=2
fsrup=4
uedown=6
ueup=8
cr=10
d0mu=0
d0=2
jpsi=4


fit=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
report='('

# [nom,  unc,  fdown,unc,  fup,  unc,  uedown,unc, ueup, unc,  CR,   unc]
for name,sample in rbList:
    #if 'crup' in name: continue
    #print sample[d0],sample[d0+1],sample[d0mu],sample[d0mu+1],sample[jpsi],sample[jpsi+1]
    #fit[rbList.index(str(name))] += sample[d0]/sample[d0+1]**2
    #fit[rbList.index(str(name))+1] += 1./sample[d0+1]**2
    #fit[rbList.index(str(name))] += sample[d0mu]/sample[d0mu+1]**2
    #fit[rbList.index(str(name))+1] += 1./sample[d0mu+1]**2
    #fit[rbList.index(str(name))] += sample[jpsi]/sample[jpsi+1]**2
    #fit[rbList.index(str(name))+1] += 1./sample[jpsi+1]**2
    fit[0] += sample[nom]/sample[nom+1]**2 
    fit[1] += 1./sample[nom+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[nom], sample[nom+1])
    idx=2
    diff = sample[nom] - sample[fsrdown]
    if diff > 0: idx=2
    else: idx=4
    fit[idx] += sample[fsrdown]/sample[fsrdown+1]**2
    fit[idx+1] += 1./sample[fsrdown+1]**2
    diff = sample[nom] - sample[fsrup]
    if diff > 0: idx=2
    else: idx=4
    fit[idx] += sample[fsrup]/sample[fsrup+1]**2
    fit[idx+1] += 1./sample[fsrup+1]**2
    diff = sample[nom] - sample[uedown]
    if diff > 0: idx=2
    else: idx=4
    fit[idx] += sample[uedown]/sample[uedown+1]**2
    fit[idx+1] += 1./sample[uedown+1]**2
    diff = sample[nom] - sample[ueup]
    if diff > 0: idx=2
    else: idx=4
    fit[idx] += sample[ueup]/sample[ueup+1]**2
    fit[idx+1] += 1./sample[ueup+1]**2
    diff = sample[nom] - sample[cr]
    if diff > 0: idx=2
    else: idx=4
    fit[idx] += sample[cr]/sample[cr+1]**2
    fit[idx+1] += 1./sample[cr+1]**2
    #fit[2] += abs(sample[fsrdown]/sample[fsrdown+1]**2)
    #fit[3] += 1./sample[fsrdown+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[fsrdown], sample[fsrdown+1])
    #fit[2] += abs(sample[uedown]/sample[uedown+1]**2)
    #fit[3] += 1./sample[uedown+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[uedown], sample[uedown+1])
    #fit[2] += abs(sample[cr]/sample[cr+1]**2)
    #fit[3] += 1./sample[cr+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[cr], sample[cr+1])
    #fit[4] += abs(sample[fsrup]/sample[fsrup+1]**2)
    #fit[5] += 1./sample[fsrup+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[fsrup], sample[fsrup+1])
    ##fit[4] += abs(sample[ueup]/sample[ueup+1]**2)
    ##fit[5] += 1./sample[ueup+1]**2 
    ##report += '%.2f / %.2f**2 + ' % (sample[ueup], sample[ueup+1])
    #fit[4] += abs(sample[cr]/sample[cr+1]**2)
    #fit[5] += 1./sample[cr+1]**2 
    #report += '%.2f / %.2f**2' % (sample[cr], sample[cr+1])
report += ') / ('
for name,sample in rbList:
    for i in xrange(0, len(sample)):
        if i%2 == 0: continue
        report += '1/%.2f**2' % sample[i]
        if i < len(sample)-1: report += ' + '
report += ')'
#print report

final=[fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3], 1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])]
print 'rB = %f +/- %f, rB_down = %f +/- %f, rB_up = %f +/- %f' % (fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3],1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5]))
print 'rB = %f +/- %f (stat) + %f (syst) - %f (syst)' % (final[0], final[1], abs(final[4]-final[0]), abs(final[2]-final[0]))
print 'rB = %.2f +/- %.2f (stat) + %.2f (syst) - %.2f (syst)' % (final[0], final[1], float(abs(final[4]-final[0])), float(abs(final[2]-final[0])))
