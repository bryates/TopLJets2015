from __future__ import print_function
import optparse
import json
import math
from collections import OrderedDict
import os

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,              type='string')

(opt, args) = parser.parse_args()

#read list of samples
jsonFile = open(opt.json,'r')
full = True if "jpT" in opt.json else False
fullstat = True

#jsonFile = open('data/era2016/rbFit.range','r')
#jsonFile = open('data/era2016/rbFit.range2','r')
#jsonFile = open('data/era2016/rbFit.pt','r')
#jsonFile = open('data/era2016/rbFit_jpT.json','r')
rbList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
jsonFile.close()


syst=['Nom', 'ISR-up', 'ISR-down', 'FSR-up', 'FSR-down', 'Underlying event up', 'Underlying event down', 'Color reconnection', 'Lepton selection up', 'Lepton selection down', 'Pile-up up', 'Pile-up down', 'Tracker efficiency up', 'Tracker efficiency down', 'Trigger selection up', 'Trigger selection down', 'JER up', 'JER down', 'JSF up', 'JSF down', '\\textsc{me}/\\textsc{ps} up', '\\textsc{me}/\\textsc{ps} down', 'Top pT', 'Top mass up', 'Top mass down']
syst=['Nom', 'ISR-up', 'ISR-down', 'FSR-up', 'FSR-down', 'Underlying event up', 'Underlying event down', 'Color reconnection', 'Lepton selection up', 'Lepton selection down', 'Pile-up up', 'Pile-up down', 'Tracker efficiency up', 'Tracker efficiency down', 'Trigger selection up', 'Trigger selection down', 'JER up', 'JER down', '\\textsc{me}/\\textsc{ps} up', '\\textsc{me}/\\textsc{ps} down', 'W up', 'W down', 'c up', 'c down', 'fit fcn']
#syst=['Nom', 'isr-up', 'isr-down', 'fsr-up', 'fsr-down', 'ueup', 'uedown', 'Color reconnetion', 'Lepton up', 'Lepton down', 'PU up', 'PU down', 'Tracker up', 'Tracker down', 'Trigger up', 'Trigger down', 'JER up', 'JER down', '\\textsc{me}/\\textsc{ps} up', '\\textsc{me}/\\textsc{ps} down']#, 'W up', 'W down', 'c up', 'c down', 'fit-fcn']
print(len(syst))
skip=['Top pT', 'Top mass up', 'Top mass down', 'JSF up', 'JSF down']#, 'JER up', 'JER down']
fit=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
sup = [[0] * len(syst)] * 3
sdown = [[0] * len(syst)] * 3
jup = [0] * len(syst)
jdown = [0] * len(syst)
dup = [0] * len(syst)
ddown = [0] * len(syst)
mup = [0] * len(syst)
mdown = [0] * len(syst)
fsr = [0] * 6 
fitunc = [0.,0.,0.]
report='('

# [nom,  unc,  fdown,unc,  fup,  unc,  uedown,unc, ueup, unc,  CR,   unc]
#i = 1
#while i < len(syst):
    #j = 2*i
    #print j
    #print '%s & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f \\\\' %(syst[(j/2)%len(syst)], rbList[d0][1][j], rbList[d0][1][j+1], rbList[d0mu][1][j], rbList[d0mu][1][j+1], rbList[jpsi][1][j], rbList[jpsi][1][j+1])
    ##print 'Nom & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f' %(rbList[0][1][i], rbList[0][1][i+1], rbList[1][1][i], rbList[1][1][i+1], rbList[2][1][i], rbList[2][1][i+1])
    #i+=1

total = 12;
for l in xrange(0,3):
    continue
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] in skip: 
            i+=1
            j+=2
            continue
        #if 'W ' in syst[i] and l is not 1: 
        #    rbList[l][1].append(rbList[l][1][0])
        #    rbList[l][1].append(rbList[l][1][1])
        #    rbList[l][1].append(rbList[l][1][0])
        #    rbList[l][1].append(rbList[l][1][1])
        #if 'c ' in syst[i] and l is not 1: 
        #    rbList[l][1].append(rbList[l][1][0])
        #    rbList[l][1].append(rbList[l][1][1])
        #    rbList[l][1].append(rbList[l][1][0])
        #    rbList[l][1].append(rbList[l][1][1])
        if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
            print('%s/down: %.4f %.4f ' % (syst[i], rbList[l][1][j], rbList[l][1][j+2]), end='\n')
            i+=1
            j+=2
        else:
            print('%s/: %.4f %.4f ' % (syst[i], rbList[l][1][j], rbList[l][1][j+2]), end='\n')
            i+=2
            j+=4
    print('\n')

print('LaTex')
print('>>>>>>')
i = 1
j = 2
up=[0.,0.,0.,0.,0.,0.]
down=[0.,0.,0.,0.,0.,0.]
#print('Nominal & ', end='')
#print('%.3f & ' % (rbList[0][1][0]), end='')
#print('%.3f & ' % (rbList[1][1][0]), end='')
#print('%.3f & \\\\\n' % (rbList[2][1][0]), end='')
#while i < len(syst):
    #if syst[i] in skip: 
      #i+=1
      #j+=2
      #continue
    #
    #if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
        #print(syst[i], ' & ', end='')
        #print('%.3f & ' % (rbList[0][1][j]), end='')
        #print('%.3f & ' % (rbList[1][1][j]), end='')
        #print('%.3f '   % (rbList[2][1][j]), end='')
        #i+=1
        #j+=2
    #else:
        #print(syst[i], ' & ', end='')
        ##print up
        #print('%.3f & ' %     (rbList[0][1][j+2]), end='')
        #print('%.3f & ' %     (rbList[1][1][j+2]), end='')
        #print('%.3f \\\\\n' % (rbList[2][1][j+2]), end='')
        #print(syst[i+1], ' & ', end='')
        ##print('%.3f$\pm$%.3f & ' % (rbList[0][1][0]-rbList[0][1][j], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        ##print down
        #print('%.3f & ' % (rbList[0][1][j]), end='')
        #print('%.3f & ' % (rbList[1][1][j]), end='')
        #print('%.3f ' %   (rbList[2][1][j]), end='')
        #i+=2
        #j+=4
    #print('\\\\')

while i < len(syst):
    if syst[i] in skip: 
        i+=1
        j+=2
        continue
    try:
        test=rbList[0][1][j]
    except:
        break
    if "Lepton" in syst[i] or "Trigger" in syst[i] or "Tracker" in syst[i]:
        fitunc[0] += (abs(rbList[0][1][0] - rbList[0][1][j]) + abs(rbList[0][1][0] - rbList[0][1][j+2]))/2.
        fitunc[1] += (abs(rbList[1][1][0] - rbList[1][1][j]) + abs(rbList[1][1][0] - rbList[1][1][j+2]))/2.
        fitunc[2] += (abs(rbList[2][1][0] - rbList[2][1][j]) + abs(rbList[2][1][0] - rbList[2][1][j+2]))/2.
    if "Top mass" in syst[i]: #FSR do not symmetrize
        #sdown[0][i] = abs(rbList[0][1][j]-rbList[0][1][0])
        #sdown[1][i] = abs(rbList[1][1][j]-rbList[1][1][0])
        #sdown[2][i] = abs(rbList[2][1][j]-rbList[2][1][0])
        down[0] = (down[0]**2 + (rbList[0][1][j]-rbList[0][1][0])**2)**.5
        down[2] = (down[2]**2 + (rbList[1][1][j]-rbList[1][1][0])**2)**.5
        down[4] = (down[4]**2 + (rbList[2][1][j]-rbList[2][1][0])**2)**.5
    if "FSR" in syst[i]: #FSR separate
        #sdown[0][i] = abs(rbList[0][1][j]-rbList[0][1][0])
        #sdown[1][i] = abs(rbList[1][1][j]-rbList[1][1][0])
        #sdown[2][i] = abs(rbList[2][1][j]-rbList[2][1][0])
        fsr[1] = rbList[0][1][j]-rbList[0][1][0]
        fsr[3] = rbList[1][1][j]-rbList[1][1][0]
        fsr[5] = rbList[2][1][j]-rbList[2][1][0]
    elif "ISR" in syst[i]: #symmetrize, ISR separate stat
        #sdown[0][i] = max(abs(rbList[0][1][j]-rbList[0][1][0]), abs(rbList[0][1][j+1]-rbList[0][1][1]))
        #sdown[1][i] = max(abs(rbList[1][1][j]-rbList[1][1][0]), abs(rbList[1][1][j+1]-rbList[1][1][1]))
        #sdown[2][i] = max(abs(rbList[2][1][j]-rbList[2][1][0]), abs(rbList[2][1][j+1]-rbList[2][1][1]))
        down[0] = (down[0]**2 + max(abs(rbList[0][1][j]-rbList[0][1][0]), abs(rbList[0][1][j+1]-rbList[0][1][1]))**2)**.5
        down[2] = (down[2]**2 + max(abs(rbList[1][1][j]-rbList[1][1][0]), abs(rbList[1][1][j+1]-rbList[1][1][1]))**2)**.5
        down[4] = (down[4]**2 + max(abs(rbList[2][1][j]-rbList[2][1][0]), abs(rbList[2][1][j+1]-rbList[2][1][1]))**2)**.5
    elif syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
        #sdown[0][i] = abs(rbList[0][1][j]-rbList[0][1][0])**2
        #sdown[1][i] = abs(rbList[1][1][j]-rbList[1][1][0])**2
        #sdown[2][i] = abs(rbList[2][1][j]-rbList[2][1][0])**2
        down[0] = (down[0]**2 + abs(rbList[0][1][j]-rbList[0][1][0])**2)**.5
        down[2] = (down[2]**2 + abs(rbList[1][1][j]-rbList[1][1][0])**2)**.5
        down[4] = (down[4]**2 + abs(rbList[2][1][j]-rbList[2][1][0])**2)**.5
    else: #symmetrize, ISR separate stat
        #sdown[0][i] = (abs(rbList[0][1][j+2]-rbList[0][1][0])**2 + abs(rbList[0][1][j]-rbList[0][1][0])**2) / 2
        #sdown[1][i] = (abs(rbList[1][1][j+2]-rbList[1][1][0])**2 + abs(rbList[1][1][j]-rbList[1][1][0])**2) / 2
        #sdown[2][i] = (abs(rbList[2][1][j+2]-rbList[2][1][0])**2 + abs(rbList[2][1][j]-rbList[2][1][0])**2) / 2
        down[0] = (down[0]**2 + (abs(rbList[0][1][j]-rbList[0][1][0])**2 + abs(rbList[0][1][j]-rbList[0][1][0])**2) / 2)**.5
        down[2] = (down[2]**2 + (abs(rbList[1][1][j]-rbList[1][1][0])**2 + abs(rbList[1][1][j]-rbList[1][1][0])**2) / 2)**.5
        down[4] = (down[4]**2 + (abs(rbList[2][1][j]-rbList[2][1][0])**2 + abs(rbList[2][1][j]-rbList[2][1][0])**2) / 2)**.5
     #down[0] = (down[0]**2 + (rbList[0][1][j]-rbList[0][1][0])**2)**.5
     #down[2] = (down[2]**2 + (rbList[1][1][j]-rbList[1][1][0])**2)**.5
     #down[4] = (down[4]**2 + (rbList[2][1][j]-rbList[2][1][0])**2)**.5
     #down[0] = (down[0]**2 + (max(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2)**.5,(rbList[0][1][j]-rbList[0][1][0]))**2))**.5
     #down[2] = (down[2]**2 + (max(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2)**.5,(rbList[1][1][j]-rbList[1][1][0]))**2))**.5
     #down[4] = (down[4]**2 + (max(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2)**.5,(rbList[2][1][j]-rbList[2][1][0]))**2))**.5
    down[1] = (down[1]**2 + abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2))**.5
    down[3] = (down[3]**2 + abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2))**.5
    down[5] = (down[5]**2 + abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2))**.5
    
    if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
        print(syst[i], ' & ', end='')
        #print('%.3f & ' %     (rbList[0][1][j]-rbList[0][1][0]), end='')
        #print('%.3f & ' %     (rbList[1][1][j]-rbList[1][1][0]), end='')
        #print('%.3f' % (rbList[2][1][j]-rbList[2][1][0]), end='')
        print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2+rbList[0][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f & ' %     (rbList[1][1][j]-rbList[1][1][0], pow(abs(rbList[1][1][j+1]**2+rbList[1][1][1]**2),0.5)), end='')
        print('%.3f$\pm$%.3f \\\\\n' % (rbList[2][1][j]-rbList[2][1][0], pow(abs(rbList[2][1][j+1]**2+rbList[2][1][1]**2),0.5)), end='')
        #print('%.3f & ' % (rbList[0][1][0]-rbList[0][1][j]), end='')
        #print('%.3f & ' % (rbList[1][1][0]-rbList[1][1][j]), end='')
        #print('%.3f '   % (rbList[2][1][0]-rbList[2][1][j]), end='')
        i+=1
        j+=2
    else:
        print(syst[i], ' & ', end='')
        #print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0], pow(abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2),0.5)), end='')
        #print up
        if "SR" in syst[i] or "Top mass" in syst[i] or full or fullstat:
            print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0], pow(abs(rbList[0][1][j+3]**2+rbList[0][1][1]**2),0.5)), end='')
            print('%.3f$\pm$%.3f & ' %     (rbList[1][1][j+2]-rbList[1][1][0], pow(abs(rbList[1][1][j+3]**2+rbList[1][1][1]**2),0.5)), end='')
            print('%.3f$\pm$%.3f \\\\\n' % (rbList[2][1][j+2]-rbList[2][1][0], pow(abs(rbList[2][1][j+3]**2+rbList[2][1][1]**2),0.5)), end='')
        else:
            print('%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0]), end='')
            print('%.3f & ' %     (rbList[1][1][j+2]-rbList[1][1][0]), end='')
            print('%.3f \\\\\n' % (rbList[2][1][j+2]-rbList[2][1][0]), end='')
        #print('%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0]), end='')
        #print('%.3f & ' %     (rbList[1][1][j+2]-rbList[1][1][0]), end='')
        #print('%.3f \\\\\n' % (rbList[2][1][j+2]-rbList[2][1][0]), end='')
        print(syst[i+1], ' & ', end='')
        #print('%.3f$\pm$%.3f & ' % (rbList[0][1][0]-rbList[0][1][j], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        #print down
        if "SR" in syst[i] or "Top mass" in syst[i] or full or fullstat:
            #print('%.3f & ' %     max(rbList[0][1][j]-rbList[0][1][0],pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
            #print('%.3f & ' %     max(rbList[1][1][j]-rbList[1][1][0],pow(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2),0.5)), end='')
            #print('%.3f ' %       max(rbList[2][1][j]-rbList[2][1][0],pow(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2),0.5)), end='')
            print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2+rbList[0][1][1]**2),0.5)), end='')
            print('%.3f$\pm$%.3f & ' %     (rbList[1][1][j]-rbList[1][1][0], pow(abs(rbList[1][1][j+1]**2+rbList[1][1][1]**2),0.5)), end='')
            print('%.3f$\pm$%.3f \\\\\n' % (rbList[2][1][j]-rbList[2][1][0], pow(abs(rbList[2][1][j+1]**2+rbList[2][1][1]**2),0.5)), end='')
        else:
            print('%.3f & ' %     (rbList[0][1][j]-rbList[0][1][0]), end='')
            print('%.3f & ' %     (rbList[1][1][j]-rbList[1][1][0]), end='')
            print('%.3f ' %       (rbList[2][1][j]-rbList[2][1][0]), end='')
        #print('%.3f$\pm$%.3f & ' % (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        #print('%.3f$\pm$%.3f & ' % (rbList[1][1][j]-rbList[1][1][0], pow(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2),0.5)), end='')
        #print('%.3f$\pm$%.3f ' %   (rbList[2][1][j]-rbList[2][1][0], pow(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2),0.5)), end='')
        #print('%.3f & ' % (rbList[0][1][j]-rbList[0][1][0]), end='')
        #print('%.3f & ' % (rbList[1][1][j]-rbList[1][1][0]), end='')
        #print('%.3f ' %   (rbList[2][1][j]-rbList[2][1][0]), end='')
        if "Top mass" in syst[i]: #FSR do not symmetrize
            sup[0][i] = abs(rbList[0][1][j+2]-rbList[0][1][0])
            sup[1][i] = abs(rbList[1][1][j+2]-rbList[1][1][0])
            sup[2][i] = abs(rbList[2][1][j+2]-rbList[2][1][0])
            up[0] = (up[0]**2 + (rbList[0][1][j+2]-rbList[0][1][0])**2)**.5
            up[2] = (up[2]**2 + (rbList[1][1][j+2]-rbList[1][1][0])**2)**.5
            up[4] = (up[4]**2 + (rbList[2][1][j+2]-rbList[2][1][0])**2)**.5
        if "FSR" in syst[i]: #FSR separate
            #sdown[0][i] = abs(rbList[0][1][j]-rbList[0][1][0])
            #sdown[1][i] = abs(rbList[1][1][j]-rbList[1][1][0])
            #sdown[2][i] = abs(rbList[2][1][j]-rbList[2][1][0])
            fsr[0] = rbList[0][1][j+2]-rbList[0][1][0]
            fsr[2] = rbList[1][1][j+2]-rbList[1][1][0]
            fsr[4] = rbList[2][1][j+2]-rbList[2][1][0]
        elif "ISR" in syst[i]: #symmetrize, ISR separate stat
            sup[0][i] = max(abs(rbList[0][1][j+2]-rbList[0][1][0]), abs(rbList[0][1][j+3]-rbList[0][1][1]))
            sup[1][i] = max(abs(rbList[1][1][j+2]-rbList[1][1][0]), abs(rbList[1][1][j+3]-rbList[1][1][1]))
            sup[2][i] = max(abs(rbList[2][1][j+2]-rbList[2][1][0]), abs(rbList[2][1][j+3]-rbList[2][1][1]))
            up[0] = (up[0]**2 + max(abs(rbList[0][1][j+2]-rbList[0][1][0]), abs(rbList[0][1][j+3]-rbList[0][1][1]))**2)**.5
            up[2] = (up[2]**2 + max(abs(rbList[1][1][j+2]-rbList[1][1][0]), abs(rbList[1][1][j+3]-rbList[1][1][1]))**2)**.5
            up[4] = (up[4]**2 + max(abs(rbList[2][1][j+2]-rbList[2][1][0]), abs(rbList[2][1][j+3]-rbList[2][1][1]))**2)**.5
        else: #symmetrize, ISR separate stat
            sup[0][i] = (abs(rbList[0][1][j+2]-rbList[0][1][0])**2 + abs(rbList[0][1][j]-rbList[0][1][0])**2) / 2
            sup[1][i] = (abs(rbList[1][1][j+2]-rbList[1][1][0])**2 + abs(rbList[1][1][j]-rbList[1][1][0])**2) / 2
            sup[2][i] = (abs(rbList[2][1][j+2]-rbList[2][1][0])**2 + abs(rbList[2][1][j]-rbList[2][1][0])**2) / 2
            up[0] = (up[0]**2 + (abs(rbList[0][1][j+2]-rbList[0][1][0])**2 + abs(rbList[0][1][j]-rbList[0][1][0])**2) / 2)**.5
            up[2] = (up[2]**2 + (abs(rbList[1][1][j+2]-rbList[1][1][0])**2 + abs(rbList[1][1][j]-rbList[1][1][0])**2) / 2)**.5
            up[4] = (up[4]**2 + (abs(rbList[2][1][j+2]-rbList[2][1][0])**2 + abs(rbList[2][1][j]-rbList[2][1][0])**2) / 2)**.5

        #up[2] = (up[2]**2 + (rbList[1][1][j+2]-rbList[1][1][0])**2)**.5
        #up[4] = (up[4]**2 + (rbList[2][1][j+2]-rbList[2][1][0])**2)**.5
        #up[0] = (up[0]**2 + (max(abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2)**.5,(rbList[0][1][j+2]-rbList[0][1][0]))**2))**.5
        #up[2] = (up[2]**2 + (max(abs(rbList[1][1][j+3]**2-rbList[1][1][1]**2)**.5,(rbList[1][1][j+2]-rbList[1][1][0]))**2))**.5
        #up[4] = (up[4]**2 + (max(abs(rbList[2][1][j+3]**2-rbList[2][1][1]**2)**.5,(rbList[2][1][j+2]-rbList[2][1][0]))**2))**.5
        up[1] = (up[1]**2 + abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2))**.5
        up[3] = (up[3]**2 + abs(rbList[1][1][j+3]**2-rbList[1][1][1]**2))**.5
        up[5] = (up[5]**2 + abs(rbList[2][1][j+3]**2-rbList[2][1][1]**2))**.5
        i+=2
        j+=4
    print('\\\\')
#print('Total up & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f$\pm$%.3f\\\\' % (up[0],up[1],up[2],up[3],up[4],up[5]))
print('Total up & %.3f & %.3f & %.3f\\\\' % (up[0],up[2],up[4]))
print('Total down & %.3f & %.3f & %.3f\\\\' % (down[0],down[2],down[4]))
print()
print('Fit procedure & %.3f & %.3f & %.3f' % (fitunc[0]/3, fitunc[1]/3, fitunc[2]/3))
print('======')
print('Don\'t forget to use sed to change 0.00 to <0.001: %s/-\{0,1\}0\.000/$<$0.001/g\n')
print('$r_{\PQb}=%0.3f \pm %0.3f \stat ^{%+0.3f}_{%+0.3f} \syst ^{%+0.3f} _{%+0.3f} \\textrm{(FSR)}$' % (rbList[0][1][0], rbList[0][1][1], up[0], -down[0], fsr[0], fsr[1]))
print('$r_{\PQb}=%0.3f \pm %0.3f \stat ^{%+0.3f}_{%+0.3f} \syst ^{%+0.3f} _{%+0.3f} \\textrm{(FSR)}$' % (rbList[1][1][0], rbList[1][1][1], up[2], -down[2], fsr[2], fsr[3]))
print('$r_{\PQb}=%0.3f \pm %0.3f \stat ^{%+0.3f}_{%+0.3f} \syst ^{%+0.3f} _{%+0.3f} \\textrm{(FSR)}$' % (rbList[2][1][0], rbList[2][1][1], up[4], -down[4], fsr[4], fsr[5]))

print('Average for BLUE (jpsi, d0, d0mu)')
blueFile = open("BLUE/rbSFCor.txt",'w')
blueFile.write("\./combine<<!\n")
blueFile.write("'2016 b-quark fragmentation'\n")
blueFile.write("-1 -1 -1 internal & minuit debug level, dependency flag\n")
blueFile.write("1 3 %d  # of observables, measurements, error classes\n" % (total))
blueFile.write("'rB'    name\n")
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'     'toppt'    'tmass'\n")
print("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'    'hdamp'     'toppt'    'tmass'")
#print("'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'")#,"'jer'")
for l in xrange(0,3):
    if l==0:
        print("'jpsi'  'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'jpsi'  'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==1:
        print("'d0'    'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0'    'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==2:
        print("'d0_mu' 'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0_mu' 'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] in skip: 
          i+=1
          j+=2
          continue
        try:
            test=rbList[0][1][j]
        except:
            break
        if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            blueFile.write('%.4f ' % low)
            i+=1
            j+=2
        else:
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            high=rbList[l][1][0]-rbList[l][1][j+2]
            highe=rbList[l][1][1]-rbList[l][1][j+3]
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % ((low+high)/2), end='')
            blueFile.write('%.4f ' % ((low+high)/2))
            i+=2
            j+=4
    print('')
    blueFile.write('\n')
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("1.0 0.0 0.0 'stat'\n")
blueFile.write("0.0 1.0 0.0\n")
blueFile.write("0.0 0.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'isr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'fsr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'ue'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'cr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'lep'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pu'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pi'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'hdamp'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'toppt'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0 'tmass'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("!\n")
blueFile.close()

print('')

print('Syst up for BLUE (jpsi, d0, d0mu)')
blueFile = open("BLUE/asyUp.txt",'w')
blueFile.write("\./combine<<!\n")
blueFile.write("'2016 b-quark fragmentation'\n")
blueFile.write("-1 -1 -1 internal & minuit debug level, dependency flag\n")
blueFile.write("1 3 %d  # of observables, measurements, error classes\n" % (total))
blueFile.write("'rB'    name\n")
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'     'toppt'    'tmass'\n")
print("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'    'hdamp'     'toppt'    'tmass'")
#print("'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'")#,"'jer'")
for l in xrange(0,3):
    if l==0:
	print("'jpsi'  'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'jpsi'  'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==1:
        print("'d0'    'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0'    'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==2:
        print("'d0_mu' 'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0_mu' 'rB' %.4f      %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] in skip: 
          i+=1
          j+=2
          continue
        try:
            test=rbList[0][1][j]
        except:
            break
        if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            blueFile.write('%.4f ' % low)
            i+=1
            j+=2
        else:
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            high=rbList[l][1][0]-rbList[l][1][j+2]
            highe=rbList[l][1][1]-rbList[l][1][j+3]
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % high, end='')
            blueFile.write('%.4f ' % high)
            i+=2
            j+=4
    print('')
    blueFile.write('\n')
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("1.0 0.0 0.0 'stat'\n")
blueFile.write("0.0 1.0 0.0\n")
blueFile.write("0.0 0.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'isr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'fsr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'ue'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'cr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'lep'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pu'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pi'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'hdamp'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'toppt'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0 'tmass'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("!\n")
blueFile.close()
print('')

print('Syst down for BLUE (jpsi, d0, d0mu)')
blueFile = open("BLUE/asyDown.txt",'w')
blueFile.write("\./combine<<!\n")
blueFile.write("'2016 b-quark fragmentation'\n")
blueFile.write("-1 -1 -1 internal & minuit debug level, dependency flag\n")
blueFile.write("1 3 %d  # of observables, measurements, error classes\n" % (total))
blueFile.write("'rB'    name\n")
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'     'toppt'    'tmass'\n")
print("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'    'hdamp'     'toppt'    'tmass'")
#print("'fsr'","'ue'","'cr'","'lep'","'trig'","'trk'","'pu'","'pi'")#,"'jer'")
for l in xrange(0,3):
    if l==0:
        print("'jpsi'  'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'jpsi'  'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==1:
        print("'d0'    'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0'    'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    if l==2:
        print("'d0_mu' 'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]), end='')
        blueFile.write("'d0_mu' 'rB' %.4f       %.4f   " % (rbList[l][1][0],rbList[l][1][1]))
    i = 1
    j = 2
    while i < len(syst):
        if syst[i] in skip: 
          i+=1
          j+=2
          continue
        try:
            test=rbList[0][1][j]
        except:
            break
        if syst[i] == 'Color reconnection' or syst[i] == 'Top pT':
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            if(low<lowe): low=lowe
            print('%.4f ' % low, end='')
            blueFile.write('%.4f ' % low)
            i+=1
            j+=2
        else:
            low=rbList[l][1][0]-rbList[l][1][j]
            lowe=rbList[l][1][1]-rbList[l][1][j+1]
            high=rbList[l][1][0]-rbList[l][1][j+2]
            highe=rbList[l][1][1]-rbList[l][1][j+3]
            if(low<lowe): low=lowe
            if(high<highe): high=highe
            print('%.4f ' % low, end='')
            blueFile.write('%.4f ' % low)
            i+=2
            j+=4
    print('')
    blueFile.write('\n')
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("\n")
blueFile.write("1.0 0.0 0.0 'stat'\n")
blueFile.write("0.0 1.0 0.0\n")
blueFile.write("0.0 0.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'isr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'fsr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'ue'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'cr'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'lep'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pu'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'pi'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'hdamp'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'toppt'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("\n")
blueFile.write("1.0 1.0 1.0 'tmass'\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("1.0 1.0 1.0\n")
blueFile.write("!\n")
blueFile.close()
#print('')

#for name,sample in rbList:
    ##if 'crup' in name: continue
    ##print sample[0],sample[0+1],sample[d0mu],sample[d0mu+1],sample[jpsi],sample[jpsi+1]
    ##fit[rbList.index(str(name))] += sample[0]/sample[0+1]**2
    ##fit[rbList.index(str(name))+1] += 1./sample[0+1]**2
    ##fit[rbList.index(str(name))] += sample[d0mu]/sample[d0mu+1]**2
    ##fit[rbList.index(str(name))+1] += 1./sample[d0mu+1]**2
    ##fit[rbList.index(str(name))] += sample[jpsi]/sample[jpsi+1]**2
    ##fit[rbList.index(str(name))+1] += 1./sample[jpsi+1]**2
    #fit[0] += sample[nom]/sample[nom+1]**2 
    #fit[1] += 1./sample[nom+1]**2 
    #report += '%.2f / %.2f**2 + ' % (sample[nom], sample[nom+1])
    #idx=2
    #diff = sample[nom] - sample[fsrdown]
    #if diff > 0: idx=2
    #else: idx=4
    #fit[idx] += sample[fsrdown]/sample[fsrdown+1]**2
    #fit[idx+1] += 1./sample[fsrdown+1]**2
    #diff = sample[nom] - sample[fsrup]
    #if diff > 0: idx=2
    #else: idx=4
    #fit[idx] += sample[fsrup]/sample[fsrup+1]**2
    #fit[idx+1] += 1./sample[fsrup+1]**2
    #diff = sample[nom] - sample[uedown]
    #if diff > 0: idx=2
    #else: idx=4
    #fit[idx] += sample[uedown]/sample[uedown+1]**2
    #fit[idx+1] += 1./sample[uedown+1]**2
    #diff = sample[nom] - sample[ueup]
    #if diff > 0: idx=2
    #else: idx=4
    #fit[idx] += sample[ueup]/sample[ueup+1]**2
    #fit[idx+1] += 1./sample[ueup+1]**2
    #diff = sample[nom] - sample[cr]
    #if diff > 0: idx=2
    #else: idx=4
    #fit[idx] += sample[cr]/sample[cr+1]**2
    #fit[idx+1] += 1./sample[cr+1]**2
    ##fit[2] += abs(sample[fsrdown]/sample[fsrdown+1]**2)
    ##fit[3] += 1./sample[fsrdown+1]**2 
    ##report += '%.2f / %.2f**2 + ' % (sample[fsrdown], sample[fsrdown+1])
    ##fit[2] += abs(sample[uedown]/sample[uedown+1]**2)
    ##fit[3] += 1./sample[uedown+1]**2 
    ##report += '%.2f / %.2f**2 + ' % (sample[uedown], sample[uedown+1])
    ##fit[2] += abs(sample[cr]/sample[cr+1]**2)
    ##fit[3] += 1./sample[cr+1]**2 
    ##report += '%.2f / %.2f**2 + ' % (sample[cr], sample[cr+1])
    ##fit[4] += abs(sample[fsrup]/sample[fsrup+1]**2)
    ##fit[5] += 1./sample[fsrup+1]**2 
    ##report += '%.2f / %.2f**2 + ' % (sample[fsrup], sample[fsrup+1])
    ###fit[4] += abs(sample[ueup]/sample[ueup+1]**2)
    ###fit[5] += 1./sample[ueup+1]**2 
    ###report += '%.2f / %.2f**2 + ' % (sample[ueup], sample[ueup+1])
    ##fit[4] += abs(sample[cr]/sample[cr+1]**2)
    ##fit[5] += 1./sample[cr+1]**2 
    ##report += '%.2f / %.2f**2' % (sample[cr], sample[cr+1])
#report += ') / ('
#for name,sample in rbList:
    #for i in xrange(0, len(sample)):
        #if i%2 == 0: continue
        #report += '1/%.2f**2' % sample[i]
        #if i < len(sample)-1: report += ' + '
#report += ')'
##print report
#
##final=[fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3], 1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])]
##print('rB = %f +/- %f, rB_down = %f +/- %f, rB_up = %f +/- %f' % (fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3],1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])))
##print('rB = %f +/- %f (stat) + %f (syst) - %f (syst)' % (final[0], final[1], abs(final[4]-final[0]), abs(final[2]-final[0])))
##print('rB = %.2f +/- %.2f (stat) + %.2f (syst) - %.2f (syst)' % (final[0], final[1], float(abs(final[4]-final[0])), float(abs(final[2]-final[0]))))
