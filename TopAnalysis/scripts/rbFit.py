import json
import math
from collections import OrderedDict

jsonFile = open('data/era2016/rbFit.json','r')
rbList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
jsonFile.close()

nom=0
fsrdown=2
fsrup=4
uedown=6
ueup=8
cr=10

fit=[0., 0., 0., 0., 0., 0.]
report='('

# [nom,  unc,  fdown,unc,  fup,  unc,  uedown,unc, ueup, unc,  CR,   unc]
for name,sample in rbList:
    fit[0] += sample[nom]/sample[nom+1]**2 
    fit[1] += 1./sample[nom+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[nom], sample[nom+1])
    fit[2] += sample[fsrdown]/sample[fsrdown+1]**2 
    fit[3] += 1./sample[fsrdown+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[fsrdown], sample[fsrdown+1])
    fit[2] += sample[uedown]/sample[uedown+1]**2 
    fit[3] += 1./sample[uedown+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[uedown], sample[uedown+1])
    fit[2] += sample[cr]/sample[cr+1]**2 
    fit[3] += 1./sample[cr+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[cr], sample[cr+1])
    fit[4] += sample[fsrup]/sample[fsrup+1]**2 
    fit[5] += 1./sample[fsrup+1]**2 
    report += '%.2f / %.2f**2 + ' % (sample[fsrup], sample[fsrup+1])
    #fit[4] += sample[ueup]/sample[ueup+1]**2 
    #fit[5] += 1./sample[ueup+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[ueup], sample[ueup+1])
    fit[4] += sample[cr]/sample[cr+1]**2 
    fit[5] += 1./sample[cr+1]**2 
    report += '%.2f / %.2f**2' % (sample[cr], sample[cr+1])
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
print 'rB = %.3f +/- %.3f (stat) + %.3f (syst) - %.3f (syst)' % (final[0], final[1], abs(final[4]-final[0]), abs(final[2]-final[0]))
print 'rB = %.2f +/- %.2f (stat) + %.3f (syst) - %.3f (syst)' % (final[0], final[1], float(abs(round(final[4],3)-round(final[0],3))), float(abs(round(final[2],3)-round(final[0],3))))
