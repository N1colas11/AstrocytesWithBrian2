
# Train MNIST
# transcribed from model-zoo/vision/mnist/

using Flux, Flux.Data.MNIST, Statistics
using Flux: onehotbatch, onecold, crossentropy
using Base.Iterators: repeated, partition
using Printf, BSON
#using CUDAapi  # should not be required if I do not have a GPU

# Load labels and images from Flux.Data.MNIST
@info("Loading data set")
train_labels = MNIST.labels()
train_imgs = MNIST.images()

# Bundle images together with labels and group in to minibatches
X_batch = Array{Float32}(undef, size(X[1])..., 1, length(idxs))
function make_minibatch(X, Y, idxs)
