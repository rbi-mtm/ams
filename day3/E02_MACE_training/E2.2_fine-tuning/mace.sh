#!/usr/bin/bash

num_threads=2
export OMP_NUM_THREADS=$num_threads
export MKL_NUM_THREADS=$num_threads

mace_run_train \
    --name="MACE" \
    --train_file="train.xyz" \
    --valid_fraction=0.05 \
    --test_file="test.xyz" \
    --E0s='average' \
    --model="MACE" \
    --hidden_irreps='8x0e' \
    --r_max=3 \
    --batch_size=4 \
    --max_num_epochs=?? \
    --swa \
    --start_swa=?? \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --energy_key="energy" \
    --forces_key="forces" \
    --restart_latest \
    --device=cpu \
    --default_dtype=float32 \
    --seed=42 \
    --foundation_model="small" \
    --multiheads_finetuning=False 

