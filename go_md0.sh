#!/bin/bash
#
# Written by Satoshi Ono
#

if [ $# -ne 1 ]; then
    echo "Usage: go_md0.sh Name"
    exit 1
fi

export PATH=$HOME/miniforge3/bin:$PATH
export BIN=$(cd $(dirname $0); pwd)
export PYTHONPATH=$BIN:$PYTHONPATH

mkdir 1
cd 1
python $BIN/make_md0.py -N $1 -n 1 -d 0
output=$(sbatch cano.sh)
echo $output
job_id=$(echo "$output"|awk '{print $4}')
cd ..

for num in `seq 2 14`
do
    mkdir $num
    cd $num
    python $BIN/make_md0.py -N $1 -n $num -d $job_id
    output=$(sbatch cano.sh)
    echo $output
    job_id=$(echo "$output"|awk '{print $4}')
   cd ..
done

cat > post_md0.sh <<EOF
#!/bin/bash -l
#SBATCH -o post.out
#SBATCH -J $1_post
#SBATCH --dependency=afterok:$job_id
export PATH=\$HOME/miniforge3/bin:\$PATH
export PYTHONPATH=$BIN:\$PYTHONPATH
python $BIN/for_first.py -i 1 2 3 4 5 6 7 8 9 10 11 12 13 14 -n md1
mv md1 ..
# find . -name md.vmcmd.ene -exec rm {} \;
EOF

sbatch post_md0.sh

cat > plot.gp <<EOF
set term dumb
plot for [i=1:14] sprintf("%d/e1.pdd",i) t sprintf("%d",i)
EOF


