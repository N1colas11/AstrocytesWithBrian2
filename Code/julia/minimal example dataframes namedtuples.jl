using DataFrames
#---------------------------------------------------------
function tst()
	# I added Pandas to the packages, and removed it, and yet, the functsions are still available.
	# Since both Pandas and DataFrames have a DataFrame method, I am forced to prepend DataFrames to
	# some or most of the DataFrames methods.
	Jleak = [0.1, 0.3]
	JIPR  = [0.2, 0.4]
	println(typeof(Jleak)) # Array{Float64,1}

	V = Vector{Any}
	types = ([V for i ∈ 1:2])

	nb_el = 2
	tuple1 = Tuple([[0.,0.] for i ∈ 1:nb_el])
	tuple2 = Tuple([V for i ∈ 1:nb_el])
	named_tuple1 = NamedTuple{(:a,:b),Tuple{Int64,Int64}}
	named_tuple2 = NamedTuple{(:a,:b),Tuple{V,V}}
	DF = DataFrames
	# The next five lines all give errors.
	#append!(DF.DataFrame(), named_tuple1)
	#push!(DF.DataFrame(), named_tuple1)  # does not work
	#append!(DF.DataFrame(), named_tuple2)
	#df = DF.DataFrame(named_tuple1)  # does not work
	#df = DF.DataFrame(named_tuple2)  # does not work

	# However, I can initiailize the DataFrame with values
	named_tuple3 = NamedTuple{(:a,:b),Tuple{1., 2.}}
	named_tuple4 = NamedTuple{(:a, :b)}([1., 2.])
	append!(DF.DataFrame(), named_tuple3)  # does not work
	append!(DF.DataFrame(), named_tuple4)
	#named_tuple4 = NamedTuple{(:a,:b),Tuple{V,V}}
end
tst()
