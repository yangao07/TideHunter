#!/usr/bin/env bash
err=(0.13 0.15 0.18 0.2)
cnt=100
period=(300 500 1000 2000 3000)
copy_num=10
var_type=NOR
seed=10

if [ $# -ne 2 ] ; then
    echo "Usage:"
    echo "$0 in.fa out_dir"
    exit
fi

in_fa=$1
out_dir=$2

mkdir $out_dir 2> /dev/null
for e in ${err[@]}
do
    for p in ${period[@]}
    do
        sim_out_dir=${out_dir}/sim_e${e}_p${p}_c${cnt}_n${copy_num}_${var_type}
        mkdir $sim_out_dir 2> /dev/null
        sim=sim.fa
        echo "python ./run_simu.py -e $e -p $p -c $cnt -n $copy_num -v $var_type -s $seed -o $sim_out_dir $in_fa $sim"
        python ./run_simu.py -e $e -p $p -c $cnt -n $copy_num -v $var_type -s $seed -o $sim_out_dir $in_fa $sim
    done
done
