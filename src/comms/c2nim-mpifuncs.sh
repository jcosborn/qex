#!/bin/sh

if ln -s mpifuncs.h mpi.h; then
  o="mpifuncs.nim"
  t="mpifuncs_.nim"
  c2nim --header -o:$t mpi.h
  rm mpi.h
  echo "# $o" >$o
  echo "# created with c2nim from mpifuncs.h" >>$o
  echo >>$o
  cat $t >>$o
  rm $t
fi
