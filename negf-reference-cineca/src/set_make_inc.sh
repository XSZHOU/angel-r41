#!/bin/bash

exitStatus=0

# for Nanohub
BASE_SRC_PATH=$1

rm -f make.inc

my_host=$HOSTNAME
if [[ "$my_host" =~ "coates" ]]; then
    echo "matched coates";
    ln -s make.inc.coates.gnu make.inc
else 
  if [[ "$my_host" =~ "nanoHUB" ]]; then
    echo "matched nanoHUB";
    sed -e "s;@SRC_PATH@;${BASE_SRC_PATH};" make.inc.nanohub > make.inc
    exitStatus=$?
  else 
    if [[ "$my_host" =~ "ncnlnx07" ]]; then
      echo "matched kubuntu";
      ln -s make.inc.kubuntu make.inc
    else
      echo "could not determine configuration. please edit set_make_inc.sh!";
      exitStatus=1
    fi
  fi
fi

exit $exitStatus
