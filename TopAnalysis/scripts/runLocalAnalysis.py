import os
import sys
import optparse
import ROOT

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',       dest='input',       help='input directory with files or single file',  default=None,       type='string')
    parser.add_option('-o', '--out',      dest='outDir',      help='output directory',                           default='analysis', type='string')
    parser.add_option(      '--ch',       dest='channel',     help='channel',                                    default=13,         type=int)
    parser.add_option(      '--charge',   dest='charge',      help='charge',                                     default=0,          type=int)
    parser.add_option(      '--flav',     dest='flav',        help='flavour splitting (for single files)',       default=0,          type=int)
    parser.add_option(      '--isTT',     dest='isTT',        help='ttbar sample (for single files)',            default=False,      action='store_true')
    parser.add_option(       '--norm',     dest='norm',        help='normalization scale (for single files)',    default=1.0,        type=float)

    (opt, args) = parser.parse_args()

    ROOT.AutoLibraryLoader.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    ROOT.gROOT.LoadMacro('src/ReadTree.cc+')
    from ROOT import ReadTree

    if '.root' in opt.input:
        ReadTree(opt.input,opt.outDir,opt.channel,opt.charge,opt.norm,opt.isTT,opt.flav)
        

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
