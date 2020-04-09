
# Train MNIST

using Flux, Flux.Data.MNIST #, Statistics
using Flux: @epochs, onehotbatch, onecold, crossentropy, throttle, mse
using Base.Iterators: repeated, partition
using Printf, BSON
#using Random
using CUDAapi
if has_cuda()
	@info "CUDA is on"
	import CuArrays
	CuArrays.allowscalar(false)
end

# Encode MNIST images ascompressed vectors that can later be decoded back into images
#

imgs = MNIST.images()

# Partition into batches of size 1000
# NOTE: Flux.data is a function, so there is a conflict
dataa = [float(hcat(vec.(imgs)...)) for imgs in partition(imgs, 1000)]
dataa = gpu.(dataa)

N = 32 # Size of the encoding (number of latent dimensions)

# You can try to make the encoder/decoder network larger.
# Also, the output of encoder is a coding of the given input.
# In this case, the input dimension is 28^2 and the output dimension
# of the encoder is 32. This implies that the coding is a compressed
# representation. We can make lossy compression via this encoder .
encoder = Dense(28^2, N, leakyrelu) |> gpu
decoder = Dense(N, 28^2, leakyrelu) |> gpu

m = Chain(encoder, decoder)

loss(x) = mse(m(x), x)

# Do not output data more than once every 5 seconds
evalcb = throttle(() -> @show(loss(dataa[1])), 5)
opt = ADAM()

# Run 10 epochs
#@epochs 10 Flux.train!(loss, params(m), zip(dataa), opt, cb = evalcb)
for epoch in 1:10
	println("==> Epoch: ", epoch)
	Flux.train!(loss, params(m), zip(dataa), opt, cb = evalcb)
end

# Sample output
#
using Images
using QuartzImageIO

img(x::Vector) = Gray.(reshape(clamp.(x, 0, 1), 28, 28))

function sample()
	# 20 random digits
	before = [imgs[i] for i in rand(1:length(imgs), 20)]
	# before and after images
	# if there no GPU, there is no need for gpu() and cpu() functions
	# vec reshapes array from 28x28 to 784
	after = img.(map(x -> cpu(m)(float(vec(x))), before))
	# Stack them all together
	hcat(vcat.(before, after)...)
end

cd(@__DIR__)
save("sample_ae.png", sample())
