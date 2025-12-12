#!/usr/bin/bash

num_threads=2
export OMP_NUM_THREADS=$num_threads
export MKL_NUM_THREADS=$num_threads

mace_run_train \
    --name="MACE" \
    --train_file="???" \
    --valid_fraction=0.05 \
    --test_file="???" \
    --E0s='average' \
    --model="MACE" \
    --num_channels=8 \
    --max_L=0 \
    --r_max=3 \
    --batch_size=4 \
    --max_num_epochs=??? \
    --swa \
    --start_swa=??? \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --energy_key="energy" \
    --forces_key="forces" \
    --restart_latest \
    --device=cpu \
    --default_dtype=float32 \
    --seed=42 

