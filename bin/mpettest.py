#!/usr/bin/env python3

import errno
import os
import os.path as osp
import pytest
import run_tests

#Run the tests and save the input args
args=run_tests.main()

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

#Create the list of arguments for pytest
pytest_args=["--baseDir="+osp.join(args.test_dir,"ref_outputs"),
             "--modDir="+args.output_dir,
             osp.join(osp.dirname(osp.abspath(__file__)),"../tests/compare_tests.py")
            ]

#Add addtional options if a test list was provided
if args.tests:
    pytest_args.append("--tests="+" ".join(args.tests))
    pytest_args.append("--skip-analytic")


#Run pytest
pytest.main(pytest_args)

#remove files if they exist 
silentremove("LA4_8rep.MWF") 
silentremove("Short_PreDiag_000173.000") 
silentremove("test_mwf_LA4.000")
