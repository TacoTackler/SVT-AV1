#!/bin/sh

CURR_DIR=$(cd $(dirname $0); pwd)
TIME=$(date "+%Y%m%d%H%M%S")
WORK_DIR=${CURR_DIR}/SVT-AV1_${TIME}
git clone https://github.com/Cidana-Developers/SVT-AV1.git ${WORK_DIR}
cd ${WORK_DIR}
git remote add upstream https://github.com/OpenVisualCloud/SVT-AV1.git
git fetch upstream
git rebase upstream/master
git checkout coverage_test
git rebase master
cd ${WORK_DIR}/code_coverage_report
time --output=coverage_test_time.txt ./generate_report.sh
cd ${WORK_DIR}/..
