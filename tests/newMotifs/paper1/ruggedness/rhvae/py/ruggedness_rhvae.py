from pythae.pipelines import TrainingPipeline
from pythae.models import RHVAE, RHVAEConfig
from pythae.trainers import BaseTrainerConfig
import numpy as np
import pandas as pd
import torch
import sys

# Run a Riemannian-Hamiltonian variational autoencoder

# Set seed
np.random.seed(42)

# Second argument is the model
model = sys.argv[1]

DATA_PATH = "/g/data/ht96/nb9894/newMotifs/paper1/ruggedness/log3/"
DATA_PATH = "/mnt/e/Documents/GitHub/SLiMTests/tests/newMotifs/paper1/ruggedness/rhvae/py/"

input_dim = 7

if (model == "FFLC1" or model == "FFLI1"):
    input_dim = 9
elif (model == "FFBH"):
    input_dim = 11


cfg = RHVAEConfig(input_dim=(input_dim,), 
                  latent_dim = 2, 
                  n_lf = 3, 
                  temperature = 0.8)

rhvae = RHVAE(cfg)

train_cfg = BaseTrainerConfig(
    output_dir = "manifold_train_" + model,
    num_epochs = 100,
    learning_rate=1e-3,
    per_device_train_batch_size=512,
    per_device_eval_batch_size=512,
    steps_saving = None
)

pipeline = TrainingPipeline(
    model = rhvae,
    training_config = train_cfg
)

# Load in dataframe
data = pd.read_csv(DATA_PATH + "d_ruggedness_" + model + ".csv")

# Load in centroids
centroids = pd.read_csv(DATA_PATH + "d_centroids_" + model + ".csv")

# Remove first column (fitness)
x = data.to_numpy()[:,1:]
c = centroids.to_numpy()[:,1:]

# Train the model
pipeline(train_data = c)

rhvae.eval()


with torch.no_grad():
    z = rhvae.cpu().encoder( # Make sure model is on cpu
        torch.tensor(x, dtype=torch.float32)
    ).embedding.numpy()

data['RH1'] = z[:, 0]
data['RH2'] = z[:, 1]
data['model'] = model

# Save output
data.to_csv(DATA_PATH + "d_ruggedness_" + model + "_rh.csv", header = True, index = False)
