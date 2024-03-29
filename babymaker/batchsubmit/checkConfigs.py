#!/bin/python

import sys, os
from glob import glob

### checks a condor submission config to see if all output files are present
### if not all present, creates a resubmission config with the missing ones
### assumes that all outputs from one config live in the same directory
### usage: python checkConfig.py <original_config>

def CheckFile(input_file):
    input_lines = []
    print 'RUNNING',input_file
    try:
        f = open(input_file,'r')
        input_lines = f.readlines()
        f.close()
    except IOError:
        print 'Could not open input file, exiting:',input_file
        exit

    found_jobs = False
    preamble = []
    input_dir = ''
    output_dir = ''
    doRebal = 0
    expected_output_files = []
    for line in input_lines:
        # save preamble at the beginning of the file, for resubmission script
        #  all lines until we find "executable="
        if not found_jobs:
            if 'arguments=' in line:
                found_jobs = True
            else:
                preamble.append(line)
                continue
        if found_jobs:
            if not 'arguments' in line:
                continue
            tokens = line.split("=")[-1].split()
            if len(input_dir) == 0:
                input_dir = tokens[1]
            if len(output_dir) == 0:
                output_dir = tokens[2]
            output_file = os.path.join(output_dir, "{}.root".format(tokens[0]))
            expected_output_files.append(output_file)

    existing_output_files = glob(output_dir+'/*.root')

    missing_files = []
    for f in expected_output_files:
        expf = f
        if not any(expf in s for s in existing_output_files):
            print 'missing file:',expf
            missing_files.append(f)

    if len(missing_files) == 0:
        print 'no missing files out of %d, hurray!' % len(expected_output_files)
        exit
    else:
        print 'missing %d/%d files, creating resubmission config' % (len(missing_files),len(expected_output_files))

        newdir = os.path.split(input_file)[0]+"_resubmit"
        os.system("mkdir -p "+newdir)
        resubmit_filename = os.path.join(newdir, os.path.basename(input_file))

        # if missing files, make a resubmission config
        outf = open(resubmit_filename,'w')
        for line in preamble:
            outf.write(line)
        for mfname in missing_files:
            mf = mfname.split("/")[-1].split(".")[0]
            outf.write('\n')
            # outf.write(executable_command)
            # outf.write('transfer_executable=True\n')
            outf.write('arguments={0} {1} {2}\n'.format(mf,input_dir,output_dir))
            outf.write('queue\n')
            outf.write('\n')

    return

def CheckDir(input_dir):
    for input_file in glob(input_dir+"/config_*.cmd"): CheckFile(input_file)
    return

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No input file or directory given.")
    else:
        args = sys.argv
        for f in args[1:]:
            if os.path.isfile(f): CheckFile(f)
            elif os.path.isdir(f): CheckDir(f)
            else: print("ERROR: {} does not exist.".format(args[1]))
