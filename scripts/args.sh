#!/bin/bash
#get arguments here to call other scripts

#calling methods required for future calls 

#make sure data isnt modified
chmod 555 ./data/*

if [ $3 == 1 ]
then
    bash ./scripts/align.sh $1 $5 $7 > ./logs/single/$1"_align.log" 2>&1 
elif [ $3 == 2 ]
then
    bash ./scripts/align_all.sh $1 $5 $7 > ./logs/pooled/$1"_align.log" 2>&1
fi

if [ $3 == 1 ]
then
    bash ./scripts/extract.sh $1 $4 $3 $6 > ./logs/single/$1"_extract.log" 2>&1
elif [ $3 == 2 ]
then
    bash ./scripts/extract.sh $1 $4 $3 $6 > ./logs/pooled/$1"_extract.log" 2>&1
fi

if [ $3 == 1 ]
then
    bash ./scripts/post_process.sh $1 $3 $8 > ./logs/single/$1"_post_process.log" 2>&1
elif [ $3 == 2 ]
then
    bash ./scripts/post_process.sh $1 $3 $8 > ./logs/pooled/$1"_post_process.log" 2>&1
fi

