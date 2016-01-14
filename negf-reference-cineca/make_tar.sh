#!/bin/bash

cd ..
rm -f angel-nanohub.tar
tar -zcvf ./angel-nanohub.tar.gz \
        ./angel-nanohub/COPYING ./angel-nanohub/COPYING.LESSER \
		./angel-nanohub/definitions.negf \
		./angel-nanohub/doc/*.[^s]* \
		./angel-nanohub/examples/*.[^s]* \
		./angel-nanohub/materials/*.[^s]* \
		./angel-nanohub/scripts/*.[^s]* \
		./angel-nanohub/Makefile.coates.gnu \
        ./angel-nanohub/*/*.cpp ./angel-nanohub/*/*.h \
		./angel-nanohub/*/*/*.cpp ./angel-nanohub/*/*/*.h \
		./angel-nanohub/*/*/*/*.cpp ./angel-nanohub/*/*/*/*.h
		
