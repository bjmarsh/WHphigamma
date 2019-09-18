import glob
import sys, os
import ROOT

ROOT.gSystem.Load("CORE/CMS3_CORE.so")
# ROOT.gROOT.ProcessLine(".L CORE/Tools/utils.cc+")
ROOT.gROOT.ProcessLine(".L WHTree/WHTree.cc+")
ROOT.gROOT.ProcessLine(".L ScanChain.C+")

tag = sys.argv[1]
indir = sys.argv[2]
infile = sys.argv[3]
samp = sys.argv[4]

if "Run2016" in indir or "Summer16" in indir:
    year = 2016
if "Run2017" in indir or "Fall17" in indir:
    year = 2017
if "Run2018" in indir or "Autumn18" in indir or "privateMC_102x" in indir:
    year = 2018

os.system("mkdir -p output_files/{0}".format(tag))

c = ROOT.TChain("Events")
for f in glob.glob("{0}/{1}.root".format(indir,infile)):
    c.Add(f)

isData = "Run20" in indir
print "isData?", isData

ROOT.ScanChain(c, tag, samp, year, isData)
