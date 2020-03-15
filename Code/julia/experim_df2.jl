using DataFrames

struct A end
struct B end

function getCurrents(::A)
           random  = rand()
           println("A")
           J_δ = [8, 6]
           J_β = [3,23]
           return (J_δ = J_δ, J_β = J_β, c = random)
end

function getCurrents(::B)
           dict = Dict()
           random  = rand()
           dict[:J_δ] = [8, 6]
           dict[:J_β] = [3,25]
           return dict
end

a = A();  output_named_tuple = getCurrents(a)
for (k,v) in zip(keys(output_named_tuple), output_named_tuple)
  # eval evaluates in the global context
  eval(:($k = $v))
end

b = B(); dict = getCurrents(b)

for (k,v) in zip(keys(dict), dict)
  # eval evaluates in the global context
  eval(:($k = $v))  # Defines J_β, J_δ
end
#------------------------
for i in 1:5
     global output_named_tuple = getCurrents()
     ac = output_named_tuple
     println(typeof(output_named_tuple))
     push!(df, output_named_tuple)
end

println(df)

for i in 1:5
     global output_named_tuple = getCurrents()
     ac = output_named_tuple
     println(typeof(output_named_tuple))
     push!(df, output_named_tuple)
end

println(df)
