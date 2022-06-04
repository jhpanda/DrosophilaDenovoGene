#! /bin/bash
## Number of nodes
#SBATCH -N 1
#SBATCH -n 4
##SBATCH -p bigmem

HOME_trRosetta="/ru-auth/local/home/jpeng/scratch/softwares/trRosetta/"
HOME_trContact="${HOME_trRosetta}/trContact"
script="${HOME_trContact}/network/predict.py"
NN_model="${HOME_trContact}/model2019_07"
fasta_dir="$2"

pref="$1"
fas="${fasta_dir}/${pref}.fasta"
inp="${pref}.a3m"
msa="${pref}_trunc.a3m"
npz="${pref}.npz"
model_num=20

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################
conda activate trRosetta

## build MSA ##
python build_MSA.py $fas

## remove headers in resulting MSA ##
python reformat.py $pref

## predict contacts ##
python $script -m $NN_model $msa $npz

## generate contact map ##
sh visual_map.sh $pref

## now deactivate trRosetta
conda deactivate

## generate models ##
mkdir -p output
cat << EOF > gen_model.sh
#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=1-${model_num}%${model_num}
ID=\${SLURM_ARRAY_TASK_ID}
pref="${pref}"
HOME_trRosetta="${HOME_trRosetta}"
tr_script="\${HOME_trRosetta}/trRosetta.py"
rsr_npz="\${pref}.npz"
fasta="\${pref}.fasta"
out="output/\${pref}_\${ID}.pdb"
python \$tr_script \$rsr_npz \$fasta \$out
EOF
#sbatch --wait gen_models.sh $pref
sbatch --wait gen_model.sh

## compute tmscore to minimum structure ##
python tmscore.py $pref
