from __future__ import print_function
import json
import math
from collections import OrderedDict
import os

#jsonFile = open('data/era2016/rbFit.range','r')
#jsonFile = open('data/era2016/rbFit.range2','r')
#jsonFile = open('data/era2016/rbFit.pt','r')
jsonFile = open('data/era2016/rbFit.json','r')
rbList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
jsonFile.close()

nom=0
isrdown=2
isrup=4
fsrdown=6
fsrup=8
uedown=8
ueup=10
cr=12
lepd=14
lepu=16
trid=18
triu=20
trkd=22
trku=24
pud=26
puu=28
jerd=30
jeru=32
pid=34
piu=36
pi=34
d0=0
d0mu=1
jpsi=2

syst=['Nom', 'ISR-up', 'ISR-down', 'FSR-up', 'FSR-down', 'Underlying event up', 'Underlying event down', 'Color reconnection', 'Lepton selection up', 'Lepton selection down', 'Pile-up up', 'Pile-up down', 'Tracker efficiency up', 'Tracker efficiency down', 'Trigger selection up', 'Trigger selection down', 'JER up', 'JER down', '\\textsc{me}/\\textsc{ps} up', '\\textsc{me}/\\textsc{ps} down']
#syst=['Nom', 'ISR-up', 'ISR-down', 'FSR-up', 'FSR-down', 'Underlying event up', 'Underlying event down', 'Color reconnection', 'Lepton selection up', 'Lepton selection down', 'Tracker efficiency up', 'Tracker efficiency down', 'Pile-up up', 'Pile-up down', '$\pi$ efficiency up', '$\pi$ efficiency down', 'Trigger selection up', 'Trigger selection down', 'JER up', 'JER down', 'hdamp up', 'hdamp down']
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

total = 9;
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
print('>>>>>>')
i = 1
j = 2
up=[0.,0.,0.,0.,0.,0.]
down=[0.,0.,0.,0.,0.,0.]
skip=['Trigger selection up', 'Trigger selection down', 'JER up', 'JER down']
while i < len(syst):
    if syst[i] in skip: 
      i+=1
      j+=2
      continue
    down[0] = (down[0]**2 + (rbList[0][1][j]-rbList[0][1][0])**2)**.5
    down[2] = (down[2]**2 + (rbList[1][1][j]-rbList[1][1][0])**2)**.5
    down[4] = (down[4]**2 + (rbList[2][1][j]-rbList[2][1][0])**2)**.5
    #down[0] = (down[0]**2 + (max(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2)**.5,(rbList[0][1][j]-rbList[0][1][0]))**2))**.5
    #down[2] = (down[2]**2 + (max(abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2)**.5,(rbList[1][1][j]-rbList[1][1][0]))**2))**.5
    #down[4] = (down[4]**2 + (max(abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2)**.5,(rbList[2][1][j]-rbList[2][1][0]))**2))**.5
    down[1] = (down[1]**2 + abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2))**.5
    down[3] = (down[3]**2 + abs(rbList[1][1][j+1]**2-rbList[1][1][1]**2))**.5
    down[5] = (down[5]**2 + abs(rbList[2][1][j+1]**2-rbList[2][1][1]**2))**.5
    
    if syst[i] == 'Color reconnection':
        print(syst[i], ' & ', end='')
        #print('%.3f$\pm$%.3f & ' % (rbList[0][1][j]-rbList[0][1][0], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        print('%.3f & ' % (rbList[0][1][0]-rbList[0][1][j]), end='')
        print('%.3f & ' % (rbList[1][1][0]-rbList[1][1][j]), end='')
        print('%.3f '   % (rbList[2][1][0]-rbList[2][1][j]), end='')
        i+=1
        j+=2
    else:
        print(syst[i], ' & ', end='')
        #print('%.3f$\pm$%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0], pow(abs(rbList[0][1][j+3]**2-rbList[0][1][1]**2),0.5)), end='')
        #print up
        print('%.3f & ' %     (rbList[0][1][j+2]-rbList[0][1][0]), end='')
        print('%.3f & ' %     (rbList[1][1][j+2]-rbList[1][1][0]), end='')
        print('%.3f \\\\\n' % (rbList[2][1][j+2]-rbList[2][1][0]), end='')
        print(syst[i+1], ' & ', end='')
        #print('%.3f$\pm$%.3f & ' % (rbList[0][1][0]-rbList[0][1][j], pow(abs(rbList[0][1][j+1]**2-rbList[0][1][1]**2),0.5)), end='')
        #print down
        print('%.3f & ' % (rbList[0][1][j]-rbList[0][1][0]), end='')
        print('%.3f & ' % (rbList[1][1][j]-rbList[1][1][0]), end='')
        print('%.3f ' %   (rbList[2][1][j]-rbList[2][1][0]), end='')
        up[0] = (up[0]**2 + (rbList[0][1][j+2]-rbList[0][1][0])**2)**.5
        up[2] = (up[2]**2 + (rbList[1][1][j+2]-rbList[1][1][0])**2)**.5
        up[4] = (up[4]**2 + (rbList[2][1][j+2]-rbList[2][1][0])**2)**.5
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
print('======')
print('Don\'t forget to use sed to change 0.00 to <0.001: %s/-\{0,1\}0\.000/$<$0.001/g\n')
print('$r_{\PQb}=%0.2f \pm %0.2f \stat ^{%0.2f}_{%0.2f} \syst$' % (rbList[0][1][0], rbList[0][1][1], max(up[0],up[1]), max(down[0],down[1])))
print('$r_{\PQb}=%0.2f \pm %0.2f \stat ^{%0.2f}_{%0.2f} \syst$' % (rbList[1][1][0], rbList[1][1][1], max(up[2],up[3]), max(down[2],down[3])))
print('$r_{\PQb}=%0.2f \pm %0.2f \stat ^{%0.2f}_{%0.2f} \syst$' % (rbList[2][1][0], rbList[2][1][1], max(up[4],up[5]), max(down[4],down[5])))

print('Average for BLUE (jpsi, d0, d0mu)')
blueFile = open("BLUE/rbSFCor.txt",'w')
blueFile.write("\./combine<<!\n")
blueFile.write("'2016 b-quark fragmentation'\n")
blueFile.write("-1 -1 -1 internal & minuit debug level, dependency flag\n")
blueFile.write("1 3 %d  # of observables, measurements, error classes\n" % (total))
blueFile.write("'rB'    name\n")
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'\n")
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
        if syst[i] == 'Color reconnection':
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
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'\n")
print("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'    'hdamp'")
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
        if syst[i] == 'Color reconnection':
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
blueFile.write("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'   'hdamp'\n")
print("                         'stat'   'isr'   'fsr'  'ue'   'cr'   'lep'  'pu'   'pi'    'hdamp'")
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
        if syst[i] == 'Color reconnection':
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
blueFile.write("!\n")
blueFile.close()
#print('')

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

#final=[fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3], 1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])]
#print('rB = %f +/- %f, rB_down = %f +/- %f, rB_up = %f +/- %f' % (fit[0]/fit[1], 1/math.sqrt(fit[1]), fit[2]/fit[3],1/math.sqrt(fit[3]), fit[4]/fit[5], 1/math.sqrt(fit[5])))
#print('rB = %f +/- %f (stat) + %f (syst) - %f (syst)' % (final[0], final[1], abs(final[4]-final[0]), abs(final[2]-final[0])))
#print('rB = %.2f +/- %.2f (stat) + %.2f (syst) - %.2f (syst)' % (final[0], final[1], float(abs(final[4]-final[0])), float(abs(final[2]-final[0]))))
