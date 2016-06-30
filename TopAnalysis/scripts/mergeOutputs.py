#! /usr/bin/env python
import os, sys
counters = {}

def isint(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

def getBaseNames(dirname):
    names = set()
    for item in os.listdir(dirname):
        filename, ext = os.path.splitext(item)
        if not ext == '.root': continue
        try:
            #if not ('MC13TeV' in filename or 'Data13TeV' in filename) : continue
            basename, number = filename.rsplit('_',1)
            print basename,number
            if not number == 'missing' and not isint(number):
                raise ValueError
            try:
                counters[basename].append(dirname+'/'+item)
            except KeyError:
                counters[basename] = [dirname+'/'+item]
            names.add(basename)

        except ValueError:
            print filename,'is single'
            names.add(filename)
    return names

try:
    inputdir = sys.argv[1]
    if not os.path.isdir(inputdir):
        print "Input directory not found:", inputdir
        exit(-1)
except IndexError:
    print "Need to provide an input directory."
    exit(-1)

noTrees=False
if len(sys.argv)>2 and sys.argv[2]=='True': noTrees=True

basenames = getBaseNames(inputdir)
print '-----------------------'
print 'Will process the following samples:', basenames

outputdir = os.path.join(inputdir)
chunkdir = os.path.join(inputdir, 'Chunks')

os.system('mkdir -p %s' % chunkdir)

for basename, files in counters.iteritems():

    filenames = " ".join(files)
    target = os.path.join(outputdir,"%s.root" % basename)

    # merging:
    print '... processing', basename
    if noTrees:
        cmd = 'hadd -f -T %s %s' % (target, filenames)
    else:
        cmd = 'hadd -f %s %s' % (target, filenames)
    os.system(cmd)

    # cleanup:
    cmd = '/bin/mv %s %s' % (filenames, chunkdir)
    os.system(cmd)


