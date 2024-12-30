# vmcmd_scripts
TTP-V-McMD scripts with patched version of gromacs-2024.4

1. Contents

add_posre.csh          Add posre in topology file
for_first.py           Make a first v-McMD run files from md0
for_next.py            Make a next v-McMD run files from previous run
for_next_P.py          Make a next v-McMD run files for production run
gen_cano.py            Generate canonical probability density function at given temparature
get_last_state.py      Get last virtual state from previous run
go_1-10P.sh            Main script to perform TTP-v-McMD runs with 10 times repeate
go_md0.sh              Script to perform pre-v-McMD run
kkkit.py               Script from omegagene (http://www.protein.osaka-u.ac.jp/rcsfp/pi/omegagene/)
make_hist.py           Make histogram file from md.vmcmd.ene
make_md0.py            Script for Prepare pre-c-McMD run
make_pdd.py            Make probability density function distributions (pdd files)
rangeinfo.py	       Script from omegagene (http://www.protein.osaka-u.ac.jp/rcsfp/pi/omegagene/)
reweight_edr.py        Reweight the trajectory to given temperature
select1000.py          Select 1000 xtc trajectories
select5000.py          Select 5000 xtc trajectories
select_xtc.py          Select n xtc trajectories

2. Demonstrate the execution method using CsA in CHCl3 as an example

# em1 Equilibration calculation
# md0 Pre-McMD run
# md1 - md9 McMD iterations
# md10 Production run

Export VMCMD_BIN directory where the scripts are placed

cd CHCl3/em1
qsub equil.sh

# Confirm that npt.gro has been created,

cd ..
mkdir md0; cd md0
$VMCMD_BIN/go_md0.sh CHCl3

# Confirm that all jobs are completed and the md1 directory has been created

$VMCMD_BIN/go_1-10P.sh >& log & # For csh

# When all jobs are completed, various files will be created in md10/All
# Among them, use 300K_fit_5000.xtc, md_Ref.gro, etc. for your analysis

3. Version information of the programs and libraries

Slurm 23.11.4

python 3.12.8
numpy  2.2.1
pandas 2.2.3
mdtraj 1.10.2
panedr 0.8.0

4. License

Scripts in bin directory are MIT License.

5. Author

Satoshi Ono (nca01750@gmail.com)
