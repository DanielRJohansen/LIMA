import os
import argparse
from typing import Any, Dict

import torch
from torch.utils.data import DataLoader, TensorDataset
import matplotlib.pyplot as plt
import itertools
import pandas as pd

import ray
from ray import tune

from Preprocessing import WaterDataset
from Model import build_encoder, build_bottleneck, build_decoder

import multiprocessing
multiprocessing.set_start_method("spawn", force=True)

# -----------------------------------------------------------------------------
# Fixed training hyperparameters
# -----------------------------------------------------------------------------
NUM_EPOCHS    = 20
BATCH_SIZE    = 128
LEARNING_RATE = 3e-4

# -----------------------------------------------------------------------------
# Composition search space (grid search)
# -----------------------------------------------------------------------------
COMPOSITION_SPACE = {
    "num_enc_blocks":   [1, 2, 3, 4],
    "enc_hidden_units": [32, 64, 128, 256],
    "enc_activation":   ["ReLU", "LeakyReLU", "ELU", "Tanh"],
    "bottleneck_type":  ["SignSTE", "GumbelBernoulli", "VQVAE"],
    "num_dec_blocks":   [1, 2, 3, 4],
    "dec_hidden_units": [32, 64, 128, 256],
    "dec_activation":   ["ReLU", "LeakyReLU", "ELU", "Tanh"]
}

COMPOSITION_SPACE_SMALL = {
    "num_enc_blocks":   [1, 2],
    "enc_hidden_units": [32],
    "enc_activation":   ["ReLU", "LeakyReLU"],
    "bottleneck_type":  ["SignSTE"],
    "num_dec_blocks":   [1, 2],
    "dec_hidden_units": [32],
    "dec_activation":   ["ReLU", "LeakyReLU"]
}

# -----------------------------------------------------------------------------
# Generator for yielding run number and val_loss
# -----------------------------------------------------------------------------
def result_generator(analysis):
    for i, trial in enumerate(analysis.trials):
        val_loss = trial.last_result.get("val_loss", None)
        yield i, val_loss


# -----------------------------------------------------------------------------
# Trainable function for Ray Tune
# -----------------------------------------------------------------------------
def train(config: Dict[str, Any], trainDataloader, validationDataloader):
    # device selection per trial
    device = torch.device("cuda")

    # build model
    encoder    = build_encoder(config).to(device)
    bottleneck = build_bottleneck(config).to(device)
    decoder    = build_decoder(config).to(device)
    model      = torch.nn.Sequential(encoder, bottleneck, decoder)

    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    loss_fn   = torch.nn.MSELoss()

    # training loop
    for _ in range(NUM_EPOCHS):
        model.train()
        for batch in trainDataloader:
            recon = model(batch)
            loss  = loss_fn(recon, batch)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    # evaluation (vectorized on GPU)
    model.eval()
    total_loss = 0.0

    errors = []
    with torch.no_grad():
        for batch in validationDataloader:
            recon = model(batch)
            total_loss += loss_fn(recon, batch).item() * batch.size(0)
            # diff shape: [B, 3 atoms, 3 coords]
            diff = (recon - batch).view(batch.size(0), 3, 3)
            # per-molecule sum of Euclidean errors over the 3 atoms
            per_mol_err = torch.norm(diff, dim=2).sum(dim=1)  # [B]
            errors.append(per_mol_err)

    # concatenate all per-molecule errors on GPU
    errors = torch.cat(errors, dim=0)  # [N_val]
    avg_loss = total_loss / len(validationDataloader.dataset)


    median_err = errors.median()
    sigma1 = errors.quantile(0.6827)
    sigma2 = errors.quantile(0.9545)
    sigma3 = errors.quantile(0.9973)

    # move only scalars to CPU/python
    median_err, err_1sigma, err_2sigma, err_3sigma = (
        median_err.item(),
        sigma1.item(),
        sigma2.item(),
        sigma3.item(),
    )

    return avg_loss, median_err, err_1sigma

def train_tune(config, train_loader, val_loader):
    avg_loss, median_err = train(config, train_loader, val_loader)
    tune.report(val_loss=avg_loss, median_error=median_err)

# -----------------------------------------------------------------------------
# Main: parallel search with live top-16 plotting
# -----------------------------------------------------------------------------

def trial_name_creator(trial):
    return f"trial_{trial.trial_id}"

if __name__ == "__main__":
    data_path = r"C:\Users\Daniel\git_repo\LIMA_data\CoordCompression\h2o\trajectory.uff"

    # Returns data on device
    train_ds, val_ds = WaterDataset(data_path, smallSet=True).Split()

    trainDataloader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    validationDataloader   = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)


    useRay = False

    if useRay:
        # initialize Ray
        ray.init(ignore_reinit_error=True)

        # convert grid to Tune config
        tune_config = {k: tune.grid_search(v) for k, v in COMPOSITION_SPACE_SMALL.items()}

        # run experiments, limit to 4 concurrent trials
        analysis = tune.run(
            tune.with_parameters(train_tune, trainDataloader=trainDataloader, validationDataloader=validationDataloader),
            config=tune_config,
            resources_per_trial={"cpu": 4, "gpu": 0.25},  # 4 trials => 1 GPU split; adjust as needed
            storage_path=r"C:\Users\Daniel\git_repo\LIMA\dev\ML\workspace",
            trial_dirname_creator=trial_name_creator,
            verbose=1,
            reuse_actors=False
        )

        # retrieve results and plot top 16
        df = analysis.results_df
    else:
        results = []
        keys = list(COMPOSITION_SPACE.keys())
        for vals in itertools.product(*COMPOSITION_SPACE.values()):
            cfg = dict(zip(keys, vals))
            avg_loss, median_err, err_1sigma = train(cfg, trainDataloader, validationDataloader)
            results.append({**cfg, "val_loss": avg_loss, "median_error": median_err})
            print("Training finished,  Val Loss: {:10.4f} Median Err {:10.4f} Sigma err: {:10.4f}".format(
                  avg_loss, median_err, err_1sigma))

        df = pd.DataFrame(results)

    top16 = df.sort_values("val_loss").head(16)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(range(16), top16["val_loss"], label="Val Loss")
    ax.plot(range(16), top16["median_error"], marker='o', linestyle='--', label="Median Error")
    ax.set_xticks(range(16))
    ax.set_xticklabels(top16.index.astype(str), rotation=45)
    ax.set_ylabel("Metric")
    ax.set_title("Top 16 Model Compositions")
    ax.legend()

    os.makedirs("./workspace/plots", exist_ok=True)
    out_path = os.path.join("./workspace/plots", "top16_results.png")
    fig.savefig(out_path)
    print(f"Saved top-16 plot to {out_path}")
