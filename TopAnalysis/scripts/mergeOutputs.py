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
            if not ('MC13TeV' in filename or 'Data13TeV' in filename) : continue
            sfilename = filename.split('_',1)
            basename, number = sfilename[1].rsplit('_',1)
            if 'filt' in filename:
                postfix=number
                basename, number=basename.rsplit('_',1)
                basename += '_'+postfix

            if not number == 'missing' and not isint(number):
                raise ValueError
            key = '_'.join([header,basename])
            try:
                counters[key].append(dirname+'/'+item)
            except KeyError:
                counters[key] = [dirname+'/'+item]
            names.add(key)

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
    #cmd = 'hadd -f %s %s' % (target, filenames)
    #os.system(cmd)

    # cleanup:
    #cmd = '/bin/mv %s %s' % (filenames, chunkdir)
    #os.system(cmd)


