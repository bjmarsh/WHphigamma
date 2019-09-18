#! /bin/bash

# # JTREE
# TAG=v3-0-0
# # for YEAR in 2017 2018; do
# for YEAR in 2016 2017 2018; do
#     for DIR in `ls -d /hadoop/cms/store/user/jguiang/rare-higgs/${YEAR}/${TAG}/*`; do
#         SAMP=`basename ${DIR}`
#         OUTDIR=/nfs-7/userdata/bemarsh/rare-higgs/${TAG}/${YEAR}
#         mkdir -p ${OUTDIR}
#         # echo "hadd ${OUTDIR}/${SAMP}.root ${DIR}/*.root"
#         # nohup nice -n19 hadd ${OUTDIR}/${SAMP}.root ${DIR}/*.root &> logs/log_${YEAR}_${SAMP}.txt &
#         echo "copyTree.py \"${DIR}/*.root\" ${OUTDIR}/${SAMP}.root -1 0 tree"
#         nohup nice -n19 copyTree.py "${DIR}/*.root" ${OUTDIR}/${SAMP}.root -1 0 tree &> logs/log_${YEAR}_${SAMP}.txt &
#     done
# done

# BTREE
TAG=v1
for YEAR in 2017; do
# for YEAR in 2016 2017 2018; do
    for DIR in `ls -d /hadoop/cms/store/user/bemarsh/WHphigamma/${TAG}/${YEAR}/*`; do
        SAMP=`basename ${DIR}`
        OUTDIR=/nfs-7/userdata/bemarsh/rare-higgs/btree/${TAG}/${YEAR}
        if [ -f ${OUTDIR}/${SAMP}.root ]; then
            continue
        fi
        mkdir -p ${OUTDIR}
        echo "copyTree.py \"${DIR}/*.root\" ${OUTDIR}/${SAMP}.root -1 0 Events \"nlep==1 && nphoton==1\""
        nohup nice -n19 copyTree.py "${DIR}/*.root" ${OUTDIR}/${SAMP}.root -1 0 Events "nlep==1 && nphoton==1" &> logs/log_${YEAR}_${SAMP}.txt &
    done
done

