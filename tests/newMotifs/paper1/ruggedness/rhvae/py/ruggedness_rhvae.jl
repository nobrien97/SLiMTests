# Code modified from Razo-Mejia et al. "Learning the shape of evolutionary landscapes:
# Geometric deep learning reveals hidden structure in phenotype-to-fitness maps"
## https://github.com/mrazomej/antibiotic_landscape
import AutoEncoderToolkit as AET
import Flux

import CSV
import DataFrames as DF

# Define model hyperparameters

# Define dimensionality of latent space
n_latent = 2
# Define number of neurons in hidden layers
n_neuron = 128

# Define RHVAE hyper-parameters
T = 0.8f0 # Temperature
lambda = 1.0f-2 # Regularization parameter
n_centroids = 256 # Number of centroids


# Define number of epochs
n_epoch = 50
# Define number of samples in batch
n_batch = 512
# Define number of samples when computing loss
n_batch_loss = 512
# Define learning rate
learning_rate = 10^-3
# Define fraction of data to be used for training
split_frac = 0.85

# Define loss function hyper-parameters
e = Float32(1E-3) # Leapfrog step size
K = 10 # Number of leapfrog steps
b = 0.3f0 # Initial temperature for tempering

# Define ELBO prefactors
logp_prefactor = [10.0f0, 0.1f0, 0.1f0]
logq_prefactor = [0.1f0, 0.1f0, 0.1f0]

# Define RHVAE hyper-parameters in a NamedTuple
rhvae_kwargs = (
    K=K,
    e=e,
    b=b,
)

# Define loss function kwargs in a NamedTuple
loss_kwargs = (
    K=K,
    e=e,
    b=b,
    logp_prefactor=logp_prefactor,
    logq_prefactor=logq_prefactor,
)

# Define by how much to subsample the time series
n_sub = 10

# Load data
for model in ["nar", "par", "fflc1", "ffli1", "ffbh"]
  # Load data
  model_name = "d_ruggedness_" .* model .* ".csv"
  d_ruggedness = CSV.read(model_name, DataFrame)
    
  # Define number of environments
  n_env = ncol(d_ruggedness) - 1

  fit_data = Matrix(d_ruggedness)
  
  # Split indexes of data into training and validation
  train_idx, val_idx = Flux.splitobs(
      1:size(fit_data, 2), at=split_frac, shuffle=true
  )
  
  # Extract train and validation data
  train_data = fit_data[:, train_idx]
  val_data = fit_data[:, val_idx]

    # Selecting centroids via k-means...

    # Select centroids via k-medoids
    centroids_data = AET.utils.centroids_kmedoids(fit_data, n_centroids)

    # Define JointGaussianLogEncoder...

    # Define encoder chain
    encoder_chain = Flux.Chain(
        # First layer
        Flux.Dense(n_env => n_neuron, Flux.identity),
        # Second layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
        # Third layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
        # Fourth layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
    )

    # Define layers for µ and log(σ)
    mu_layer = Flux.Dense(n_neuron => n_latent, Flux.identity)
    log_sigma_layer = Flux.Dense(n_neuron => n_latent, Flux.identity)

    # build encoder
    encoder = AET.JointGaussianLogEncoder(encoder_chain, mu_layer, log_sigma_layer)

    # Define SimpleGaussianDecoder...

    # Initialize decoder
    decoder = AET.SimpleGaussianDecoder(
        Flux.Chain(
            # First layer
            Flux.Dense(n_latent => n_neuron, Flux.identity),
            # Second Layer
            Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
            # Third layer
            Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
            # Fourth layer
            Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
            # Output layer
            Flux.Dense(n_neuron => n_env, Flux.identity)
        )
    )

    # Define MetricChain (learns Riemannian metric tensor of latent space)

    # Define mlp chain
    mlp_chain = Flux.Chain(
        # First layer
        Flux.Dense(n_env => n_neuron, Flux.identity),
        # Second layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
        # Third layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
        # Fourth layer
        Flux.Dense(n_neuron => n_neuron, Flux.leakyrelu),
    )

    # Define layers for the diagonal and lower triangular part of the covariance
    # matrix
    diag = Flux.Dense(n_neuron => n_latent, Flux.identity)
    lower = Flux.Dense(
        n_neuron => n_latent * (n_latent - 1) / 2, Flux.identity
    )

    # Build metric chain
    metric_chain = AET.RHVAEs.MetricChain(mlp_chain, diag, lower)

    # Define RHVAE model
    rhvae = AET.RHVAEs.RHVAE(
    encoder * decoder,
    metric_chain,
    centroids_data,
    T,
    lambda
    )

println("Checking previous model states...")

# List previous model parameters
model_states = sort(Glob.glob("$(state_dir)/beta-rhvae_epoch*.jld2"[2:end], "/"))

# Check if model states exist
if length(model_states) > 0
    # Load model state
    model_state = JLD2.load(model_states[end])["model_state"]
    # Input parameters to model
    Flux.loadmodel!(rhvae, model_state)
    # Update metric parameters
    AET.RHVAEs.update_metric!(rhvae)
    # Extract epoch number
    epoch_init = parse(
        Int, match(r"epoch(\d+)", model_states[end]).captures[1]
    ) + 1
else
    epoch_init = 1
end # if

println("Initial epoch: $epoch_init")

println("Uploading model to GPU...")

# Check if CUDA is available
if CUDA.functional()
    # Upload model to GPU
    rhvae = Flux.gpu(rhvae)
    # Upload data to GPU
    train_data = Flux.gpu(train_data)
    val_data = Flux.gpu(val_data)
end

# Explicit setup of optimizer
opt_rhvae = Flux.Train.setup(
    Flux.Optimisers.Adam(η),
    rhvae
)

println("\nTraining RHVAE...\n")

# Loop through number of epochs
for epoch in epoch_init:n_epoch
    # Define number of batches
    num_batches = size(train_data, 2) / n_batch
    # Shuffle data indexes
    idx_shuffle = Random.shuffle(1:size(train_data, 2))
    # Split indexes into batches
    idx_batches = IterTools.partition(idx_shuffle, n_batch)
    # Loop through batches
    for (i, idx_tuple) in enumerate(idx_batches)
        println("Epoch: $(epoch) | Batch: $(i) / $(length(idx_batches))")
        # Extract indexes
        idx_batch = collect(idx_tuple)
        # Train RHVAE
        loss_epoch = AET.RHVAEs.train!(
            rhvae, train_data[:, idx_batch], opt_rhvae;
            loss_kwargs=loss_kwargs, verbose=false, loss_return=true
        )
        println("Loss: $(loss_epoch)")
    end # for train_loader

    # Sample train data
    train_sample = train_data[
        :,
        StatsBase.sample(1:size(train_data, 2), n_batch_loss, replace=false)
    ]
    # Sample val data
    val_sample = val_data

    println("Computing loss in training and validation data...")
    loss_train = AET.RHVAEs.loss(rhvae, train_sample; loss_kwargs...)
    loss_val = AET.RHVAEs.loss(rhvae, val_sample; loss_kwargs...)

    # Forward pass sample through model
    println("Computing MSE in training and validation data...")
    out_train = rhvae(train_sample; rhvae_kwargs...).μ
    mse_train = Flux.mse(train_sample, out_train)
    out_val = rhvae(val_sample; rhvae_kwargs...).μ
    mse_val = Flux.mse(val_sample, out_val)

    println(
        "\n Epoch: $(epoch) / $(n_epoch)\n " *
        "   - loss_train: $(loss_train)\n" *
        "   - loss_val: $(loss_val)\n" *
        "   - mse_train: $(mse_train)\n" *
        "   - mse_val: $(mse_val)\n"
    )

    # Save checkpoint
    JLD2.jldsave(
        "$(state_dir)/beta-rhvae_epoch$(lpad(epoch, 5, "0")).jld2",
        model_state=Flux.state(rhvae) |> Flux.cpu,
        loss_train=loss_train,
        loss_val=loss_val,
        mse_train=mse_train,
        mse_val=mse_val,
        train_idx=train_idx,
        val_idx=val_idx,
    )
end # for n_epoch

# Save training data to dataframe for plotting
# Loading trained model...

# Find model file
model_file = first(Glob.glob("$(out_dir)/rhvae_model*.jld2"))
# List epoch parameters
model_states = Glob.glob("$(state_dir)/*.jld2")

# Initialize dataframe to store files metadata
df_meta = DF.DataFrame()

# Loop over files
for f in model_states
    # Extract epoch number from file name using regular expression
    epoch = parse(Int, match(r"epoch(\d+)", f).captures[1])
    # Load model_state file
    f_load = JLD2.load(f)
    # Extract values
    loss_train = f_load["loss_train"]
    loss_val = f_load["loss_val"]
    mse_train = f_load["mse_train"]
    mse_val = f_load["mse_val"]
    # Generate temporary dataframe to store metadata
    df_tmp = DF.DataFrame(
        :epoch => epoch,
        :loss_train => loss_train,
        :loss_val => loss_val,
        :mse_train => mse_train,
        :mse_val => mse_val,
        :model_file => model_file,
        :model_state => f,
    )
    # Append temporary dataframe to main dataframe
    global df_meta = DF.vcat(df_meta, df_tmp)
end # for f in model_states

CSV.Save(df_meta, "d_rhvae_" .* model .* ".csv")

end
## =============================================================================


# Now fit fitness profiles to the latent space

println("Load model...\n")

# Load model
rhvae = JLD2.load("$(vae_dir)/rhvae_model.jld2")["model"]
# Load parameters
model_state = JLD2.load("$(vae_dir)/rhvae_model.jld2")["model_state"]
# Input parameters to model
Flux.loadmodel!(rhvae, model_state)
# Update metric parameters
AET.RHVAEs.update_metric!(rhvae)


# Define latent space dimensions
latent = DD.Dim{:latent}([:latent1, :latent2])

# Map data to latent space
dd_latent = DD.DimArray(
    dropdims(
        mapslices(slice -> rhvae.vae.encoder(slice).μ,
            log_fitnotype_std.data,
            dims=[5]);
        dims=1
    ),
    (log_fitnotype_std.dims[2:4]..., latent, log_fitnotype_std.dims[6]),
)