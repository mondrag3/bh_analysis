 #!/usr/bin/env python

import getopt
import sys
import re
from array import array
import os
import math
import ROOT


def process(params):
    # we want to handle Ctrl+C
    sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
    sh.Add()
    sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

    # create a chain
    chain = ROOT.TChain('t3')
    for name in params.inputs:
        chain.Add(name)

    out = None

    n = 0
    weights = 0
    weights2 = 0
    for event in chain:
        name = "%d %d -> %s" % (event.id1, event.id2, ' '.join([str(x) for x in event.kf]))
        #print name, event.weight, event.me_wgt
        weights += event.weight
        weights2 += event.weight**2
        n += 1

    print "XS ", n, weights/n, math.sqrt(weights2 - weights**2/n)/(n-1)


def usage():
    print """\
Usage: analyze [OPTION...] [FILE]
Reweight events
  -o, --output='name.root'  Output name
  -d, --debug               debug

Other options:
  -h, --help                show this help message
"""


class Params:
    def __init__(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "o:dh",
                                 ["output=", "debug", "help"])
        except getopt.GetoptError, err:
            print str(err)
            usage()
            sys.exit(2)

        self.output = 'test.root'
        self.debug = False

        for op, oparg in opts:
            if op in ("-h", "--help"):
                usage()
                sys.exit()
            elif op in ("-o", "--output"):
                self.output = oparg
            elif op in ("-d", "--debug"):
                self.debug = True
            else:
                assert False, "unhandled option"

        if len(args) > 0:
            self.inputs = args[:]
        else:
            print "Error: missing input files"
            usage()
            sys.exit(2)


def main():
    params = Params()
    process(params)


if __name__ == '__main__':
    main()
