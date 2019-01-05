#!/bin/bash
if [ $# -ne 8 ] ; then
    echo "Usage:"
    echo "$0 in.fa out_prefix type(CLR/CCS) depth len-min len-max acc-min acc-max"
    exit
fi

pbsim=pbsim
cur_dir=$(pwd)
m_clr=$cur_dir/model_qc_clr
m_ccs=$cur_dir/model_qc_ccs

in_fa=$1
out_pre=$2
t=$3
if [ "$t" == "CLR" ]; then
    m=$m_clr
else
    m=$m_ccs
fi
depth=$4
min_l=$5
max_l=$6
min_acc=$7
max_acc=$8

echo "$pbsim --prefix $out_pre --data-type $t --length-min $min_l --length-mean $min_l --depth $depth --length-max $max_l --accuracy-min $min_acc --accuracy-max $max_acc --model_qc $m $in_fa"
$pbsim --prefix $out_pre --data-type $t --length-min $min_l --length-mean $min_l --depth $depth --length-max $max_l --accuracy-min $min_acc --accuracy-max $max_acc --model_qc $m $in_fa


