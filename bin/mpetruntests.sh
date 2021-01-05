#!/bin/bash
set -x
set -e
rm -rf workdir
mkdir workdir
cd workdir

mkdir stable_branch
if [ "$1" = "stable" ]; then
  git --work-tree=stable_branch/ checkout $1 -- mpet
else
  git --work-tree=stable_branch/ checkout remotes/origin/$1 -- mpet
fi

cp ../run_tests.py .
ln -s ../../mpet .
ln -s ../../tests .

#run tests for current commit
coverage run --source=../../mpet/ run_tests.py --test_dir ./tests --output_dir ../../bin/workdir/modified > /dev/null
COVERALLS_REPO_TOKEN=$COVERALLS_REPO_TOKEN coveralls || : #Dont fret if if fails

#run tests for stable branch
ln -sf stable_branch/mpet .
python run_tests.py --test_dir ./tests --output_dir ../../bin/workdir/stable > /dev/null

