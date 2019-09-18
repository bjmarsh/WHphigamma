#! /bin/bash

tar -cJf input.tar.xz ../do.py ../ScanChain.C ../vhAngles.h ../WHTree/*.cc ../WHTree/*.h ../CORE/*.h ../CORE/Tools/*.h ../CORE/Tools/*.cc \
    ../CORE/*.so ../jetCorrections ../jsons ../CORE/Tools/datasetinfo/*.h ../CORE/Tools/jetcorr/*.h \
    -C /home/users/bemarsh/scripts/ copyTree.py
