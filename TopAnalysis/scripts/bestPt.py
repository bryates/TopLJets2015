import json
import math
from collections import OrderedDict
import os

def rnd(num, base):
    return base * round(float(num)/base)

jsonFile = open('data/era2016/rbFit.json','r')
rbList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
jsonFile.close()
systs=['172v5', 'isr-down', 'isr-up', 'fsr-down', 'fsr-up', 'uedown', 'ueup', 'erdON','down_LEP', 'up_LEP', 'down_PU', 'up_PU', 'down_PI', 'up_PI', 'down_TRIGGER', 'up_TRIGGER', 'down_JER', 'up_JER', 'hdampdown', 'hdampup', 'tpt']

print rbList
for r in range(len(systs)):
    rb = rbList[1][1][r]
    if r > 0: rb = rbList[0][1][r+r%2]
    syst = systs[r]
    if syst == 'erdON': syst = 'GluonMove_erdON'
    rbr = int(rnd(rb,0.025)*1000)
    if rbr%50==0 and rbr%100!=0: rbr += 5
    cmd = "cp LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/mcVdata_" + syst + "_" + str(rbr) + "_BCDEFGH_d0.png LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/best/mcVdata_" + syst + "_best_d0.png"
    print cmd
    os.system(cmd)
    #print rnd(float(rb), 0.025)
for r in range(len(systs)):
    rb = rbList[2][1][r]
    if r > 0: rb = rbList[1][1][r+r%2]
    syst = systs[r]
    rbr = int(rnd(rb,0.025)*1000)
    if rbr%50==0 and rbr%100!=0: rbr += 5
    print rb,rbr
    cmd = "cp LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/mcVdata_" + syst + "_" + str(rbr) + "_d0mu_tag_BCDEFGH_d0mu.png LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/best/mcVdata_" + syst + "_best_d0mu_tag_d0mu.png"
    print cmd
    os.system(cmd)
    #print rnd(float(rb), 0.025)

for r in range(len(systs)):
    rb = rbList[0][1][r]
    if r > 0: rb = rbList[2][1][r+r%2]
    syst = systs[r]
    rbr = int(rnd(rb,0.025)*1000)
    if rbr%50==0 and rbr%100!=0: rbr += 5
    cmd = "cp LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/mcVdata_" + syst + "_" + str(rbr) + "_BCDEFGH_jpsi.png LJets2015/2016/mtop/www/meson/morph/jpT/mcVdata/best/mcVdata_" + syst + "_best_jpsi.png"
    print cmd
    os.system(cmd)
    #print rnd(float(rb), 0.025)
