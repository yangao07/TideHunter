## Simulation datasets used in [TideHunter paper](https://doi.org/10.1093/bioinformatics/btz376)

15 simulation datasets were included in `simu.tar.gz`. 
```
tar -zxvf simu.tar.gz
```
There are three folders, each one corresponds to one type of the evaluation performed in the paper:
```
copy_num/
err_rate/
repeat_size/
```

In each of the three folders, 5 different copy numbers/error rates/repeat sizes were used to generate 1000 raw long reads with tandem repeats: `sim.fa`. The folder name indicates the simulation parameters:
```
copy_num/sim_e0.15_s1000_c10/sim.fa
copy_num/sim_e0.15_s1000_c20/sim.fa
copy_num/sim_e0.15_s1000_c2/sim.fa
copy_num/sim_e0.15_s1000_c3/sim.fa
copy_num/sim_e0.15_s1000_c5/sim.fa
err_rate/sim_e0.13_s1000_c10/sim.fa
err_rate/sim_e0.15a_s1000_c10/sim.fa
err_rate/sim_e0.15b_s1000_c10/sim.fa
err_rate/sim_e0.16_s1000_c10/sim.fa
err_rate/sim_e0.20_s1000_c10/sim.fa
repeat_size/sim_e0.15_s1000_c10/sim.fa
repeat_size/sim_e0.15_s100_c10/sim.fa
repeat_size/sim_e0.15_s2000_c10/sim.fa
repeat_size/sim_e0.15_s3000_c10/sim.fa
repeat_size/sim_e0.15_s500_c10/sim.fa

```
where,
```
e: error rate
s: repeat size
c: copy number
```

In each folder, in addition to the simulated long reads `sim.fa`, a concatemer of two copies of the repeat unit sequence was also generated: `sim.fa.tr`. The two-copy concatemer can be used to calculate the accuracy of the called consensus sequence.
```
copy_num/sim_e0.15_s1000_c10/sim.fa.tr
copy_num/sim_e0.15_s1000_c20/sim.fa.tr
copy_num/sim_e0.15_s1000_c2/sim.fa.tr
copy_num/sim_e0.15_s1000_c3/sim.fa.tr
copy_num/sim_e0.15_s1000_c5/sim.fa.tr
err_rate/sim_e0.13_s1000_c10/sim.fa.tr
err_rate/sim_e0.15a_s1000_c10/sim.fa.tr
err_rate/sim_e0.15b_s1000_c10/sim.fa.tr
err_rate/sim_e0.16_s1000_c10/sim.fa.tr
err_rate/sim_e0.20_s1000_c10/sim.fa.tr
repeat_size/sim_e0.15_s1000_c10/sim.fa.tr
repeat_size/sim_e0.15_s100_c10/sim.fa.tr
repeat_size/sim_e0.15_s2000_c10/sim.fa.tr
repeat_size/sim_e0.15_s3000_c10/sim.fa.tr
repeat_size/sim_e0.15_s500_c10/sim.fa.tr
```
