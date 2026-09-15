from pythae.pipelines import TrainingPipeline
from pythae.models import RHVAE, RHVAEConfig
from pythae.trainers import BaseTrainerConfig
import numpy as np
import pandas as pd
import torch
import sys

# Run a Riemannian-Hamiltonian variational autoencoder

# Second argument is the model
model = sys.argv[1]

DATA_PATH = "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/log3/"

cfg = RHVAEConfig(input_dim=(7,), 
                  latent_dim = 2, 
                  n_lf = 3, 
                  temperature = 0.8)

rhvae = RHVAE(cfg)

train_cfg = BaseTrainerConfig(
    output_dir = 'test_model_nar',
    num_epochs = 50,
    learning_rate=1e-3,
    per_device_train_batch_size=128,
    per_device_eval_batch_size=128,
    steps_saving = None
)

pipeline = TrainingPipeline(
    model = rhvae,
    training_config = train_cfg
)

# Load in dataframe
data = pd.read_csv(DATA_PATH + "d_ruggedness_" + model + ".csv")

# Remove first column (fitness)
x = data.to_numpy()[:,1:]

# Train on a smaller subset of x
x_sbst = x[np.random.choice(x.shape[0], 30000, replace = False), :]

pipeline(train_data = x_sbst)

rhvae.eval()

with torch.no_grad():
    z = rhvae.encoder(
        torch.tensor(x, dtype=torch.float32)
    ).embedding.cpu().numpy()

data['RH1'] = z[:, 0]
data['RH2'] = z[:, 1]
data['model'] = model

# Save output
data.to_csv(DATA_PATH + "d_ruggedness_" + model + "_rh.csv", header = False, index = False)
