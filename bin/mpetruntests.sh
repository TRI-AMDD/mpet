#!/bin/bash
set -x
set -e
rm -rf workdir
mkdir workdir
cd workdir

mkdir stable_branch
cp -r ../../.git stable_branch/
cd stable_branch
git branch --all
git reset remotes/origin/$1 --hard
cd ../

cp ../run_tests.py .
ln -s ../../mpet .
ln -s ../../tests .

#run tests for current commit
coverage run --source=../../mpet/ run_tests.py ../../bin/workdir/modified > /dev/null
COVERALLS_REPO_TOKEN=$COVERALLS_REPO_TOKEN coveralls || : #Dont fret if if fails

#run tests for stable branch
ln -sf stable_branch/mpet .
python run_tests.py ../../bin/workdir/stable > /dev/null

