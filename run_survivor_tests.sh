source venv/bin/activate
source ssshtest
# Override truv to not use coverage if it's missing, but it should be in venv
# Actually setup_test.sh sets truv="coverage run ..."
# Let's just make sure coverage is there.
TESTSRC=repo_utils/
source $TESTSRC/setup_test.sh
# If coverage is not installed, we can override truv
if ! command -v coverage &> /dev/null; then
    truv="python3 -m truvari.__main__"
fi
rm -rf $OD
mkdir -p $OD
source $TESTSRC/sub_tests/survivor.sh
