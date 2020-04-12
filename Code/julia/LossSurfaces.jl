using Flux
using BSON: @save, @load
using CuArrays

batchsize = 64

dat_x = Flux.Data.FashionMNIST.images(:train)
dat_x = [reshape(Float32.(x), 28, 28, 1, 1) for x in dat_x]
dat_x = gpu.([reduce((x, y)->cat(x, y, dims=4), batch)
            for batch in Iterators.partition(dat_x, batchsize)])

dat_y = Flux.Data.FashionMNIST.labels(:train)
dat_y = gpu.([Flux.onehotbatch(batch, 0:9)
            for batch in Iterators.partition(dat_y, batchsize)])

dat = zip(dat_x, dat_y);
#------------------------------------------------
function myTrain!(loss, model, data, opt)
    ps = params(model)
    for d in data
        gs = gradient(ps) do
            training_loss = loss(d...)
            println("loss: ", training_loss)
            return training_loss
        end
    end
end

#-------------------------------------------------------
function large_vgg()
    vgg = cpu(Chain(
        Conv((3,3), 1  =>64,  relu, pad=(1,1)), BatchNorm(64), # 28x28
        Conv((3,3), 64 =>64,  relu, pad=(1,1)), BatchNorm(64), # 28x28
        # Change maxpool function signature (upgraded Flux)
        MaxPool((2,2)),  # 28x28 -> 14x14
        Conv((3,3), 64 =>128, relu, pad=(1,1)), BatchNorm(128),
        Conv((3,3), 128=>128, relu, pad=(1,1)), BatchNorm(128),
        MaxPool((2,2)), # 14x14 -> 7x7
        Conv((3,3), 128=>256, relu, pad=(1,1)), BatchNorm(256),
        Conv((3,3), 256=>256, relu, pad=(1,1)), BatchNorm(256),
        Conv((3,3), 256=>256, relu, pad=(1,1)), BatchNorm(256),
        MaxPool((2,2)), # 7x7 -> 3x3
        Conv((3,3), 256=>512, relu, pad=(1,1)), BatchNorm(512),
        Conv((3,3), 512=>512, relu, pad=(1,1)), BatchNorm(512),
        Conv((3,3), 512=>512, relu, pad=(1,1)), BatchNorm(512),
        MaxPool((2,2)),  # 3x3 -> 1x1
        x -> reshape(x, :, size(x,4)),
        Dense(512,  4096, relu),
        Dense(4096, 4096, relu), Dropout(0.5),
        Dense(4096, 10),
        softmax))
    @save "large_vgg" vgg
end

@info "before large_vgg()"
#large_vgg()

function run_large_vgg()
    # Loading must occur from the CPU
    @load "large_vgg" vgg
    @info "loaded large_vgg"
    vgg = gpu(vgg)
    println("After @load")
    loss(x, y) = Flux.crossentropy(vgg(x), y)
    println("Before myTrain")
    myTrain!(loss, vgg, dat, ADAM())
end

@info "before run_large_vgg()"
run_large_vgg()

#-------------------------------------------------
function small_vgg(batch_size)
    n = 32
    vgg = cpu(Chain(
        Conv((3,3), 1 => 64,  relu, pad=(1,1)), BatchNorm(64),  # 28x28=784
        Conv((3,3), 1 => 64,  relu, pad=(1,1)), BatchNorm(64),  # 28x28=784
        # size(x,4) is the batch size. : is all the other dimensions
        MaxPool((2,2)), x -> reshape(x, :, size(x,4)), # 14x14=196
        Dense(196*batch_size,  n, relu),  # 14x14
        Dense(n, n, relu), Dropout(0.5),
        Dense(n, 10),
        softmax))
        @show(vgg)
    @save "small_vgg" vgg
end
######################################
######################################
# dat is global
function run_small_vgg()
    @load "small_vgg" vgg
    loss(x, y) = begin
        println("sizes: ", size(x), size(y))
        Flux.crossentropy(vgg(x), y)
    end
    myTrain!(loss, vgg, dat, ADAM())
end

small_vgg(batchsize)
run_small_vgg()
#----------------------------------------------------------------------
