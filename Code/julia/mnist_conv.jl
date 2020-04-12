
# Train MNIST
# transcribed from model-zoo/vision/mnist/

using Flux, Flux.Data.MNIST, Statistics
using Flux: onehotbatch, onecold, crossentropy, throttle
using Base.Iterators: repeated, partition
using Printf, BSON
using Random
using CUDAapi
#using CUDAapi  # should not be required if I do not have a GPU

# Load labels and images from Flux.Data.MNIST
@info("Loading data set")
train_labels = MNIST.labels()
train_imgs = MNIST.images()

# Bundle images together with labels and group in to minibatches
function make_minibatch(X, Y, idxs)
        X_batch = Array{Float32}(undef, size(X[1])..., 1, length(idxs))
        for i in 1:length(idxs)
                X_batch[:,:,:,i] = Float32.(X[idxs[i]])
        end
        Y_batch = onehotbatch(Y[idxs], 0:9)
        return (X_batch, Y_batch)
end

batch_size = 128
mb_idxs = partition(1:length(train_imgs), batch_size)
train_set = [make_minibatch(train_imgs, train_labels, i) for i in mb_idxs]

#Prepare test set as one giant minibatch:
test_imgs = MNIST.images(:test)
test_labels = MNIST.labels(:test)

# Need the same shuffle
perm_index = shuffle(1:length(test_labels))
# GE: Make sure that test_imgs and test_labels are permuted the same way
test_imgs, test_labels = test_imgs[perm_index], test_labels[perm_index]

# GE: try and understand the float. (vectorized). Why doesn't @. work?
tX = cat(float.(test_imgs[1:1000])..., dims=4) |> gpu
tY = onehotbatch((test_labels[1:1000]), 0:9) |> gpu
test_set = make_minibatch(test_imgs, test_labels, 1:length(test_imgs))

# Define model. We use a simple convolutional architecture with
# three iterations of conv -> ReLU -> MaxPool, followed by a final Dense
# layer that feeds into a softmax probability output
@info("Constructing model...")
model = Chain(
        # First convolution, operating upon a 28x28 image
        Conv((3,3), 1=>16, pad=(1,1), relu),
        MaxPool((2,2)),

        # Second convolution, operating upon a 14x14 image
        Conv((3,3), 16=>32, pad=(1,1), relu),
        MaxPool((2,2)),

        # Third convolution, operating upon a 7x7 image
        Conv((3,3), 32=>32, pad=(1,1), relu),
        MaxPool((2,2)),

        # Reshape 3d tensor into a 2D one, at this point it should be (3,3,32,N)
        # which is where we get the 288 in the 'Dense' layer below:

        x -> reshape(x,:,size(x,4)),
        Dense(288, 10),

        # Finaly, softmax to get nice probabilities
        softmax,
)

# Load model and datasets onto GPU, if enabled
train_set = gpu.(train_set)
test_set = gpu.(test_set)
model = gpu(model)

# Make sure our model is nicely precompiled before starting our training loop
model(train_set[1][1])

# We augment `x` a little bit here, adding in random noise
augment(x) = x .+ gpu.(0.1f0 .* randn(eltype(x), size(x))) # no error
#augment(x) = @. x + gpu(0.1f0 * randn(eltype(x), size(x))) # error

paramvec(m) = vcat(map(p->reshape(p, :), params(m))...)
anynan(x) = any(isnan.(x))

# `loss()` calculates the crossentropy loss between our prediction `y_hat`
#(calculated from `model(x)`) and the gorund truth `y` . We augment the data
# a bit, adding ggaussian random noise to our image to make it more robust.
function loss(x, y)
       xhat = augment(x)
       yhat = model(xhat)
       return crossentropy(yhat, y)
end
accuracy(x, y) = mean(onecold(cpu(model(x))) .== onecold(cpu(y)))

# Train our model with the given training set using the ADAM optimizer and
# printing out performance against the test set as we go
opt = ADAM(0.001)

# callback for training: do not call callback more than once every 10 seconds
#evalcb = throttle(() -> @show(accuracy(tX, tY)), 10)
# No constraints on the frequency callback, so callback every batch
evalcb = () -> @show(accuracy(tX, tY))

@info("Beginning training loop...")

best_acc = 0.0
last_improvement = 0
for epoch_idx = 1:100
        println("==================================")
        println("==> Epoch ", epoch_idx)  # epoch marker
        global best_acc, last_improvement

        #Flux.train!(loss, params(model), train_set, opt, cb=evalcb)
		## First argument to throttle should be a function with no arguments
        Flux.train!(loss, params(model), train_set, opt, cb=throttle(evalcb, 10))

        println("before anynan")
        if anynan(paramvec(model))
                @error "NaN params"
                break
        end

        # Calculate acuracy
        println("before accuracy")
        acc = accuracy(test_set...)
        @info(@sprintf("[%d]: Test accuracy: %.4f", epoch_idx, acc))

        # If our accuracy is good enough, quit out
        if acc >= 0.999
                @info("->Early exiting: We reached our target accuracy of 99.9%")
                break
        end

        # If this is the best accuracy seen so far, save the model output
        println("before accuracy check")
        if acc >= best_acc
                @info(" -> New best accuracy! Saving model out to mnist_conv.bson")
                BSON.@save joinpath(dirname(@__FILE__), "mnist_conv.bson") params=cpu.(params(model)) epoch_idx acc
                best_acc = acc
                last_improvement = epoch_idx
        end

        # If we haven't seen improvment in 5 epochs, drop our learning rate:
        if epoch_idx - last_improvement >= 5 && opt.eta > 1.e-6
                opt.eta /= 10.0
                @warn(" -> Haven't improved in a while, dropping learning rate to $(opt.eta)!")

                # After dropping learning rate, give it a few epochs to improve
                last_improvement = epoch_idx
        end

        if epoch_idx - last_improvement >= 10
                @warn(" -> We're calling this converged.")
                break
        end
end
