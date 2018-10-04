#!/bin/sh
# sbatch --job-name=MakeMeth --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./SPLIT_TRAIN_VALI.sh
python -u split_train_vali.py