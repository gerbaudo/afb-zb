```
1103  2013-10-30 13:56:15 cd madgraph/
1106  2013-10-30 13:57:51 diff /gdata/atlas/cshimmin/gen_monohiggs/MadGraph5/madgraph/interface/madevent_interface.py interface/madevent_interface.py
1107  2013-10-30 13:58:32 cp /gdata/atlas/cshimmin/gen_monohiggs/MadGraph5/madgraph/interface/madevent_interface.py interface/madevent_interface.py
1108  2013-10-30 13:59:12 diff /gdata/atlas/cshimmin/gen_monohiggs/MadGraph5/madgraph/various/cluster.py various/cluster.py
1109  2013-10-30 13:59:22 cd ..
1113  2013-10-30 13:59:36 mkdir ../templates
1114  2013-10-30 13:59:44 ./bin/mg5

cd MadGraph5_v1_5_12
./bin/mg5
mg5>generate p p > Z , Z > l+ l-
mg5>output ../templates/test_zll



1115  2013-10-30 14:00:06 cd ..
1117  2013-10-30 14:00:07 cd templates/
1118  2013-10-30 14:00:07 ls
1119  2013-10-30 14:00:08 cd test_zz4l/
1120  2013-10-30 14:00:09 ls
1121  2013-10-30 14:00:12 vim Cards/
1122  2013-10-30 14:01:09 ls
1123  2013-10-30 14:01:32 vim Cards/me5_configuration.txt
1125  2013-10-30 14:02:36 ssh cshimmin@compute-9-28 "ls /scratch"
1129  2013-10-30 14:04:09 ls SubProcesses/
1130  2013-10-30 14:04:12 ls SubProcesses/P0_qq_zz_z_ll_z_ll/
1133  2013-10-30 14:05:13 qstat -q atlas
1134  2013-10-30 14:05:26 ./bin/generate_events --cluster
1135  2013-10-30 14:08:50 ls
1136  2013-10-30 14:08:55 ls
1137  2013-10-30 14:08:56 ls EV
1138  2013-10-30 14:09:00 ls Events/
1139  2013-10-30 14:09:09 links crossx.html
1140  2013-10-30 14:09:27 ./bin/madevent
1141  2013-10-30 14:46:42 pwd
1142  2013-10-30 14:47:37 ls
1143  2013-10-30 14:47:38 pwd
1144  2013-10-30 14:47:39 ls ..
1145  2013-10-30 14:47:42 ls ../..
1146  2013-10-30 14:47:51 pwd
1147  2013-10-30 14:47:53 ls ../
1148  2013-10-30 14:48:09 history > hist.txt



cd templates
cd test_zz4l
vim Cards/me5_configuration.txt
# see mods
# find /gdata/atlas/gerbaudo/md5/ -name me5_configuration.txt -exec grep cluster_type {} \; -print
# cluster_type = pbs
# /gdata/atlas/gerbaudo/md5/templates/test_zz4l/Cards/me5_configuration.txt
# #cluster_type = condor

test_zz4l $ ./bin/madevent
MGME5>pythia run_01 --cluster
Which programs do you want to run?
    0 / auto    : running existing card
    1 / pythia  : Pythia
    2 / pgs     : Pythia + PGS
    3 / delphes  : Pythia + Delphes.
>3

MGME5>pythia run_01 --cluster

links test_zll/crossx.html
```

Note to self
------------
By default madgraph does not consider the `b` and `b~` within `p`.
One needs to enable it by modifying

```
$ diff input/multiparticles_default.txt input/multiparticles_default.txt~
15,16c15,16
< p = g u c d s b u~ c~ d~ s~ b~
< j = g u c d s b u~ c~ d~ s~ b~
---
> p = g u c d s u~ c~ d~ s~
> j = g u c d s u~ c~ d~ s~
```

Generate sm_mirror
------------------

```
cd MadGraph5_v1_5_12
./bin/mg5
mg5> import sm_mirror
mg5> generate p p > b mu mu~
mg5> output ../templates/bm_zb_mm


mg5>import sm_mirror
mg5>generate p p > b~ mu+ mu-
mg5>output ../templates/bm_zbbar_mm

cd ../
 cp -rp templates/bm_zb_mm runs/bm_zb_mm_take1
 cp -rp templates/bm_zbbar_mm runs/bm_zbbar_mm_take1


cd runs/bm_zb_mm_take1
patch -p1 -R Cards/me5_configuration.txt ../me5_configuration.txt.patch

cd runs/bm_zbbar_mm_take1
patch -p1 -R Cards/me5_configuration.txt ../me5_configuration.txt.patch

cp  ../sm_zb_mm_take1/run_gen_phy_del.sh .
./run_gen_phy_del.sh
```
