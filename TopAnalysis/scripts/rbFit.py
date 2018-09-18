from __future__ import print_function
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
lepd=12
lepu=14
trid=16
triu=18
trkd=20
trku=22
pud=24
puu=26
jerd=28
jeru=30
pid=32
piu=34
pi=32
d0=0
d0mu=1
jpsi=2

syst=['Nom', 'FSR-up', 'FSR-down', 'Underlying event up', 'Underlying event down', 'Color reconnection', 'Lepton selection up', 'Lepton selection down', 'Trigger selection up', 'Trigger selection down', 'Tracker efficiency up', 'Tracker efficiency down', 'Pile-up up', 'Pile-up down', '$\pi$ efficiency up', '$\pi$ efficiency down', 'JER up', 'JER down']
fit=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
report='('

# [nom,  unc,  fdown,unc,  fup,  unc,  uedown,unc, ueup, unc,  CR,   unc]
#i = 1
#while i < len(syst):
    #j = 2*i
    #print j
    #print '%s & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f \\\\' %(syst[(j/2)%len(syst)], rbList[d0][1][j], rbList[d0][1][j+1], rbList[d0mu][1][j], rbList[d0mu][1][j+1], rbList[jpsi][1][j], rbList[jpsi][1][j+1])
    ##print 'Nom & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f' %(rbList[0][1][i], rbList[0][1][i+1], rbList[1][1][i], rbList[1][1][i+1], rbList[2][1][i], rbList[2][1][i+1])
    #i+=1

for l in xrange(0,3):
    continue
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] == 'Color reconnection':
            print('%.4f %.4f ' % (rbList[l][1][j], rbList[l][1][j+2]), end='')
            i+=1
            j+=2
        else:
            print('%.4f %.4f ' % (rbList[l][1][j], rbList[l][1][j+2]), end='')
            i+=2
            j+=4
    print('\n')

print('LaTex')
i = 1
j = 2
sys=[0.,0.,0.,0.,0.,0.]
while i < len(syst)-1:
    sys[0] = (sys[0]**2 + (rbList[0][1][j]-rbList[0][1][0])**2)**.5
    sys[2] = (sys[2]**2 + (rbList[1][1][j]-rbList[1][1][0])**2)**.5
    sys[4] = (sys[4]**2 + (rbList[2][1][j]-rbList[2][1][0])**2)**.5
    sys[1] = (sys[1]**2 + abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2))**.5
    sys[3] = (sys[3]**2 + abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2))**.5
    sys[5] = (sys[5]**2 + abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2))**.5
    
    if syst[i] == 'Color reconnection':
        print(syst[i], ' & ', end='')
        print('%.3f$\pm$%.3f & ' % (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f & ' % (rbList[1][1][j]-rbList[1][1][0], pow(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f '   % (rbList[2][1][j]-rbList[2][1][0], pow(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2),0.5)), end='')
        i+=1
        j+=2
    else:
        print(syst[i], ' & ', end='')
        print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0], pow(abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f & ' %     (rbList[1][1][j+2]-rbList[1][1][0], pow(abs(rbList[1][1][j+3]**2-rbList[1][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f \\\\\n' % (rbList[2][1][j+2]-rbList[2][1][0], pow(abs(rbList[2][1][j+3]**2-rbList[2][1][1]**2),0.5)), end='')
        print(syst[i+1], ' & ', end='')
        print('%.3f$\pm$%.3f & ' % (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f & ' % (rbList[1][1][j]-rbList[1][1][0], pow(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f ' %   (rbList[2][1][j]-rbList[2][1][0], pow(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2),0.5)), end='')
        sys[0] = (sys[0]**2 + (rbList[0][1][j+2]-rbList[0][1][0])**2)**.5
        sys[2] = (sys[2]**2 + (rbList[1][1][j+2]-rbList[1][1][0])**2)**.5
        sys[4] = (sys[4]**2 + (rbList[2][1][j+2]-rbList[2][1][0])**2)**.5
        sys[1] = (sys[1]**2 + abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2))**.5
        sys[3] = (sys[3]**2 + abs(rbList[1][1][j+3]**2-rbList[1][1][1]**2))**.5
        sys[5] = (sys[5]**2 + abs(rbList[2][1][j+3]**2-rbList[2][1][1]**2))**.5
        i+=2
        j+=4
    print('\\\\')
print('Total & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f' % (sys[0],sys[1],sys[2],sys[3],sys[4],sys[5]))

print('Average for BLUE (jpsi, d0, d0mu)')
print("'stat'","'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'","'jer'")
for l in xrange(0,3):
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] == 'Color reconnection':
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            i+=1
            j+=2
        else:
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            high=abs(rbList[l][1][0]-abs(rbList[l][1][j+2]))
            highe=abs(rbList[l][1][1]-abs(rbList[l][1][j+3]))
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % ((low+high)/2), end='')
            i+=2
            j+=4
    print('')

print('Syst up')
print('Average for BLUE (jpsi, d0, d0mu)')
print("'stat'","'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'","'jer'")
for l in xrange(0,3):
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] == 'Color reconnection':
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            i+=1
            j+=2
        else:
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            high=abs(rbList[l][1][0]-abs(rbList[l][1][j+2]))
            highe=abs(rbList[l][1][1]-abs(rbList[l][1][j+3]))
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % high, end='')
            i+=2
            j+=4
    print('')

print('Syst down')
print('Average for BLUE (jpsi, d0, d0mu)')
print("'stat'","'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'","'jer'")
for l in xrange(0,3):
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] == 'Color reconnection':
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            i+=1
            j+=2
        else:
            low=abs(rbList[l][1][0]-abs(rbList[l][1][j]))
            lowe=abs(rbList[l][1][1]-abs(rbList[l][1][j+1]))
            high=abs(rbList[l][1][0]-abs(rbList[l][1][j+2]))
            highe=abs(rbList[l][1][1]-abs(rbList[l][1][j+3]))
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % low, end='')
            i+=2
            j+=4
    print('')

for name,sample in rbList:
    #if 'crup' in name: continue
    #print sample[0],sample[0+1],sample[d0mu],sample[d0mu+1],sample[jpsi],sample[jpsi+1]
    #fit[rbList.index(str(name))] += sample[0]/sample[0+1]**2
    #fit[rbList.index(str(name))+1] += 1./sample[0+1]**2
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
print('rB = %f +/- %f, rB_down = %f +/- %f, rB_up = %f +/- %f' % (fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3],1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])))
print('rB = %f +/- %f (stat) + %f (syst) - %f (syst)' % (final[0], final[1], abs(final[4]-final[0]), abs(final[2]-final[0])))
print('rB = %.2f +/- %.2f (stat) + %.2f (syst) - %.2f (syst)' % (final[0], final[1], float(abs(final[4]-final[0])), float(abs(final[2]-final[0]))))
