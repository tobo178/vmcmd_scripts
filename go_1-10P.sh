#!/bin/bash
#
# Written by Satoshi Ono
#

if [ $# -ne 1 ]; then
    echo "Usage: go_1-10P.sh Name"
    exit 1
fi

export PATH=$HOME/miniforge3/bin:$PATH
export BIN=$(cd $(dirname $0); pwd)
export PYTHONPATH=$BIN:$PYTHONPATH

source $HOME/gmx2024.4_vmcmd_cpu/bin/GMXRC.bash

echo `pwd`

wait_for_completion() {
  local job_id=$1
  echo "Waiting for job $job_id to complete..."
  while squeue -j $job_id > /dev/null 2>&1 && squeue -j $job_id | grep -q "$job_id"; do
    sleep 30
  done
  echo "Job $job_id completed."
}

# 1
if [ -d md1 ]; then
  cd md1
  echo 'md1'
  job_id=$(sbatch -J "$1_md1" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md1 -n md2
fi

# 2
if [ -d md2 ]; then
  cd md2
  echo 'md2'
  job_id=$(sbatch -J "$1_md2" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md2 -n md3
fi

# 3
if [ -d md3 ]; then
  cd md3
  echo 'md3'
  sed -i 's/^20000$/10000/' temp_s
  job_id=$(sbatch -J "$1_md3" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md3 -n md4
fi

# 4
if [ -d md4 ]; then
  cd md4
  echo 'md4'
  job_id=$(sbatch -J "$1_md4" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md4 -n md5
fi

# 5
if [ -d md5 ]; then
  cd md5
  echo 'md5'
  sed -i 's/^10000$/5000/' temp_s
  job_id=$(sbatch -J "$1_md5" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md5 -n md6
fi

# 6
if [ -d md6 ]; then
  cd md6
  echo 'md6'
  job_id=$(sbatch -J "$1_md6" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md6 -n md7
fi

# 7
if [ -d md7 ]; then
  cd md7
  echo 'md7'
  job_id=$(sbatch -J "$1_md7" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md7 -n md8
fi

# 8
if [ -d md8 ]; then
  cd md8
  echo 'md8'
  job_id=$(sbatch -J "$1_md8" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next.py -c md8 -n md9
fi

# 9
if [ -d md9 ]; then
  cd md9
  echo 'md9'
  job_id=$(sbatch -J "$1_md9" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  cd ..
  python $BIN/for_next_P.py -c md9 -n md10
fi

# 10P

if [ -d md10 ]; then
  cd md10
  echo 'md10'
  job_id=$(sbatch -J "$1_md10" ./multi.sh | awk '{print $NF}')
  wait_for_completion $job_id
  python $BIN/make_pdd.py
  (cd ..;python $BIN/gen_cano.py -c md10 -o gencano.md10 -T 300)
  mkdir All
  gmx trjcat -cat -o All/md_noSolpbc.xtc -f [1-9]/md_noSolpbc.xtc [1-9][0-9]/md_noSolpbc.xtc [1-3][0-9][0-9]/md_noSolpbc.xtc > /dev/null 2>&1
  gmx eneconv -o All/md_all.edr -f [1-9]/md_nt.edr [1-9][0-9]/md_nt.edr [1-3][0-9][0-9]/md_nt.edr > /dev/null 2>&1

  cd All
  cp ../../gencano.md10/P_E_T_300.dat .
  cp ../1/md_Ref.gro .
  /usr/bin/echo -e "1\n1" | gmx trjconv -quiet -s md_Ref.gro -f md_noSolpbc.xtc -o md_fit_all.xtc -fit "rot+trans"

  echo "Potential"|gmx energy -f md_all.edr -o potential.xvg > /dev/null 2>&1
  cat potential.xvg |sed -e "/^#/d" |sed -e "/^@/d"|awk '{printf("%13.5f\n",$2)}'> compd.eto

  python $BIN/reweight_edr.py

  echo -e "1\n1" | gmx trjconv -s md_Ref.gro -f 300K.xtc -o 300K_fit.xtc -fit "rot+trans" > /dev/null 2>&1

  python $BIN/select5000.py

  cd ..
fi
