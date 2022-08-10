### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 5dcb2929-115e-459c-b98d-43ae7bcabd3a
using Pkg; Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/jclouds")

# ╔═╡ a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
begin
using Markdown
using HDF5
using DataFrames
using StatsBase
using Plots
using StatsBase
using LinearAlgebra
using Statistics
using PlutoUI
using Random
using Base
using Distributions
end

# ╔═╡ 1751666d-d26f-4d8e-b5d6-b36854dd7462
#Pkg.add("Distributions")

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


Test NB for clouds in Julia


J.A. Hernado,

Santiago, August 2022

---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly()

# ╔═╡ ceaff28e-ea37-4a92-a812-210fd2b91fad


# ╔═╡ 94bb1cff-3f57-4af6-b177-ff732cc78429


# ╔═╡ 40a8cad8-0a7b-41c5-9807-a5143c30a76b
begin
function _random()
	size = 1000000
	bins = 10
	d    = Normal()
	xbins = -5.:1:5.
	ybins = 0:0.25:1.
	ys   = rand(size) 
	xs   = rand(d, size)
	h    = fit(Histogram, (xs, ys), (xbins, ybins))
end
end

# ╔═╡ e4e0ac6b-ad55-4a96-81e7-c5ec8df53be6
begin
nn = 10
zz = 0:nn-1
zz = vcat(zz, reverse(zz))
zz = ones(nn)' .* zz
println(Base.size(zz))
xx = 0.5 .+ 0:2*nn
yy = 0.5 .+ 0:nn
xm = ones(nn)' .* xx
ym = yy'       .* ones(2*nn)
xbins = 0:1:2*nn
ybins = 0:1:nn
end

# ╔═╡ 6faf5d14-f3f2-4549-bc1b-7d46771eb1da
zz

# ╔═╡ b05e2ea0-66c7-4c40-a998-92d7b8baa37f
ybins

# ╔═╡ b7fb1cb1-941f-4bfa-ad96-7addfd1e60c0
hh = fit(Histogram, (vec(xm), vec(ym)), weights(vec(zz)), (xbins, ybins))

# ╔═╡ 5cc3372f-f9df-4063-9e4c-dfa4e43f316c
plot(hh)

# ╔═╡ 8f49762d-9e98-4fed-b55c-17d9425d0fff
histogram2d(vec(xm), vec(ym), weights = vec(zz), nbins = (xbins, ybins))

# ╔═╡ b9095b27-8f5d-40c5-9860-b7895a9955fb
"""
Converts a histogram2d into a NamedTuple with x, y, weights
"""
function histo_to_data_(h::Histogram)::NamedTuple
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in h.edges]
	xmesh   = centers[1]' .* ones(length(centers[2]))
	ymesh   = ones(length(centers[1]))' .* centers[2]
	data    = (x = vec(xmesh), y = vec(ymesh), weights = copy(vec(h.weights')))
end

# ╔═╡ 602dc3e2-bd48-4322-8687-42a03e1d0dd7
function histo_to_data(h::Histogram)
	nx, ny  = Base.size(h.weights)
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in h.edges]
	xm = ones(ny)'   .* centers[1]
	ym = centers[2]' .* ones(nx)
	zz = deepcopy(h.weights)
	return (x = vec(xm), y = vec(ym), weights = vec(zz))
end


# ╔═╡ b06f35cd-368a-4eec-b21d-68f70afe59fa
data = histo_to_data(hh)

# ╔═╡ df7b2cd3-5b2f-4b36-a200-aaedb38bd678


# ╔═╡ 62867032-c903-4b5b-b6db-51c02ddc6e6d
histogram2d(data.x, data.y, weights = data.weights, bins = (xbins, ybins))

# ╔═╡ 33108e40-21a0-4ab7-bbc3-7827cb9c7a45
md"""

## Creating delta histograms

"""

# ╔═╡ c6b36bd0-97a1-4c02-8dcf-3d1f208299ba
hh.edges[2]

# ╔═╡ 6f354604-3e53-4a6f-956b-eac04f01c949
begin
	moves  = [[i, j] for i in -1:1:1 for j in -1:1:1]
	imove0 = [i for (i, move) in enumerate(moves) if move == [0, 0]][1]
end

# ╔═╡ 1d414c02-b275-4c27-a18f-99518f02442d
hh.edges

# ╔═╡ 50e73c5e-bc43-4cf6-b085-e1890e706e4a
begin
function _deltas(h::Histogram)
	xstep, ystep = (edge[2] - edge[1] for edge in h.edges) 
	xbins, ybins = h.edges
	data = histo_to_data(h)
	xstep, ystep = (edge[2] - edge[1] for edge in h.edges) 
	println(xstep, ", ", ystep)
	his    = [fit(Histogram, (data.x .- mx * xstep, data.y .- my * ystep), weights(data.weights), (xbins, ybins)) for (mx, my) in moves]
	h0     = deepcopy(his[imove0])
	deltas = [h.weights .- h0.weights for h in his]
	dsteps  = [norm([xstep, ystep] .* move) for move in moves]
	dsteps[imove0] = 1.
	deltas = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]
	return his, deltas
end
end

# ╔═╡ ce15bb7e-e9a1-4589-97b5-c209e4525f37
his, deltas = _deltas(hh)

# ╔═╡ 98eabc22-3748-4b71-ae6b-44e98a6be83d
plot(his, layout = length(his))

# ╔═╡ 211410f2-f516-4538-aa86-78aeac2d3b19
plot(deltas, layout = length(deltas))

# ╔═╡ 5b0a5a5c-099d-45f4-ac0c-1c24bff169d0
begin
function _grad(deltas)
	nx, ny = Base.size(deltas[imove0])
	d0     = deltas[imove0]
	grad  = deepcopy(d0)
	igrad = imove0 .* ones(Int, nx, ny)
	for (i, di) in enumerate(deltas)
		#mod    = norm(moves[i])
		#dip           = di ./mod
		imask         = di .> grad
		grad[imask]  .= di[imask]
		igrad[imask] .= i
	end
	return grad, igrad
end
end

# ╔═╡ b5415946-1e85-432f-8a2e-c72eec552ee2
norm(moves[1])

# ╔═╡ 6fc37709-d5ff-475f-877a-0b48bc54f02f
grad, igrad = _grad(deltas)

# ╔═╡ f8f35d6b-fa66-4cd4-96ea-7e260fe2633e
begin
hgrad, higrad = deepcopy(hh), deepcopy(hh)
hgrad.weights = grad; higrad.weights = igrad
end

# ╔═╡ 59c3d73d-2841-4340-abb0-4accfeb0bba9
plot(higrad)

# ╔═╡ 1bf32c73-58db-46d2-b3d0-a6f2032ffafd
minimum(igrad)

# ╔═╡ 69b2e4fb-7a76-4ddd-a172-0f2677b08648
begin
xstep, ystep = 1., 1.
dmoves = Dict(1:9 .=> moves); dmoves[0] = [0, 0]
uus = [dmoves[k] for k in vec(igrad)]
us = [xi * 0.8 * xstep for (xi, yi) in uus]
vs = [yi * 0.8 * ystep for (xi, yi) in uus]	
end

# ╔═╡ 8c574ffb-dd95-4bb4-b0a3-e24d371b3c32
dmoves[8], dmoves[5], dmoves[2]

# ╔═╡ e204b639-cb6b-48d3-861b-83278bcd98b9
uus

# ╔═╡ 465ccc2b-5d58-4a9a-9a56-24be0ceb971a
vs

# ╔═╡ d9c5642d-a93c-4094-b9b3-44b821ec20c8
begin
quiver(data.x, data.y, quiver = (us, vs))
plot!(higrad, alpha = 0.3)
end

# ╔═╡ f524364b-1dd5-4445-b96d-b8f4c98504fd
md"""
Curvature
"""

# ╔═╡ 9da43e9a-d1cd-4ad5-9876-c41023f46fa0
dot(moves[1], moves[2])

# ╔═╡ b5d78c84-1584-426e-9afa-c7d15297bedf
begin
	omoves = Dict()
	imoves = deepcopy(moves)
	imoves = filter(move -> move != imoves[imove0], moves)
	for i in 1:9
		vals = 
		omoves[moves[i]] = [imove for imove in imoves 
			if dot(imove, moves[i]) == 0]
		
	end
end

# ╔═╡ d21157e9-7da5-449c-bb3a-de76084d0f62
omoves[[0, 0]]

# ╔═╡ b049b708-934f-4531-89b8-66173c47d08c
begin
	ddeltas = Dict()
	for i in 1:9
		ddeltas[moves[i]] = deltas[i]
	end
end

# ╔═╡ f8b2e785-50f4-4af2-b63c-85eeedba31ed
ddeltas

# ╔═╡ 6d5f5671-2c6f-460a-bb3b-4b1788f5efe4
omoves

# ╔═╡ 2339e4b0-e95b-4026-be77-e4306db471a5
begin
	dels = Dict()
	for imove in moves
		dels[imove] = reduce(.+, [ddeltas[move] for move in omoves[imove]])
	end
end

# ╔═╡ 3930b27b-c1e8-4886-acbc-fc3e0099c5cb
dels

# ╔═╡ e9f7ff07-54c5-4f3a-a1d6-25d4d54499a2
begin
imove = [0, 0]
heatmap(dels[imove]')
end

# ╔═╡ 01cdd761-d68a-4563-914e-97e85b6438cb
begin
function _curvature(dels)
	nx, ny = Base.size(dels[[0, 0]])
	curmin  =  1e6 .* ones(nx, ny)
	icurmin = imove0 .* ones(Int, nx, ny) 
	curmax  = -1e6 .* ones(nx, ny)
	icurmax = imove0 .* ones(Int, nx, ny) 
	for (i, move) in enumerate(moves)
		dd = dels[move]
		mask1 = dd .> curmax
		curmax[mask1] .= dd[mask1]
		icurmax[mask1] .= i
		mask2 = dd .< curmin
		curmin[mask2] .= dd[mask2]
		icurmin[mask2] .= i 
	end
	return curmax, icurmax, curmin, icurmin
end
end

# ╔═╡ cfb9266e-cfc0-45ef-af6f-1ed5be7ce9f6
curmax, icurmax, curmin, icurmin = _curvature(dels)

# ╔═╡ 8d50c759-262c-4e98-aa22-7e597a545e86


# ╔═╡ f248c208-912a-40cb-9da5-c757af5301cc
begin
function _mesh(edges)
	nx, ny = length(edges[1]), length(edges[2])
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]
	xm = ones(ny)'   .* centers[1]
	ym = centers[2]' .* ones(nx)
	return vec(xm), vec(ym)
end
	
function _quiver(idir)
	uus = [dmoves[k] for k in vec(idir)]
	us = [xi * 0.8 * xstep for (xi, yi) in uus]
	vs = [yi * 0.8 * ystep for (xi, yi) in uus]	
	return us, vs
end
end

# ╔═╡ 7f16716e-e452-4db8-8c6c-e600a0b50944


# ╔═╡ 2df93823-002f-46f7-9d64-b625a66fa688
keys(dels)

# ╔═╡ c22e56f9-1bab-499e-a80d-58afccc2ef8f
begin
xxm, yym = _mesh(hh.edges)
wws, vvs   = _quiver(icurmin)   
end

# ╔═╡ a2f8b0b4-83ba-4c41-9a54-6f56d28febde
plot(heatmap(curmax'), heatmap(icurmax'))

# ╔═╡ eebf5cc7-58d0-44b4-9c2d-b65587bbff88
plot(heatmap(curmin'), heatmap(icurmin'))

# ╔═╡ fa808f07-96b6-4e1e-83e1-0694c969823c
omoves[[1, 1]], omoves[[-1, -1]]

# ╔═╡ f9f106c6-14ff-46d0-b060-0a552c3a082c
begin
	imove_ = [0, 0]
	heatmap(dels[imove_]')
end

# ╔═╡ e12fb4a4-f24e-43db-9d12-deebc6f2b0ba
moves[4]

# ╔═╡ c005b1d2-4ecf-4b44-bfcb-549ea746707d


# ╔═╡ 3b5730e5-145c-4fed-99b0-2e1da1982f68
moves[1], moves[2]

# ╔═╡ 8f34afb7-f955-46d7-81ce-0eb282f23993
quiver

# ╔═╡ ba0ba2a7-aa5b-4026-adde-3f767baafbe1
higrad.weights

# ╔═╡ b5c836bb-0155-41bd-88a7-70e62f2fcd3a
moves[1]

# ╔═╡ 5276e76f-5af5-49c3-80b7-b177f4eab65b
minimum(higrad.weights)

# ╔═╡ eb9b5214-cfb9-4812-9e64-5b94ea419ca9
[dmoves[i] for i in 0:9]

# ╔═╡ a36ef6d9-b23d-4896-b930-11b026e762dc
begin
quiver(data.x, data.y, quiver = (us, vs))
plot!(higrad, alpha = 0.3)
end

# ╔═╡ 8c9d468a-ebf7-4873-b358-0aec1511cbe2


# ╔═╡ 0a7f02b0-e8b1-41fa-8bd9-741a0a1bcf11
plot(his[6])

# ╔═╡ d174e2cd-d37f-48ac-8960-7dc2f815ddb5
begin
igdata  = histo_to_data(higrad)
#histogram(vec(igdata.weights), bins = 10)
#vals = [dmoves[i] for i in igdata.weights]
end

# ╔═╡ f91f8c7a-63e3-4518-bc8d-7bbd8c113c5b
histogram2d(igdata.x, igdata.y, weights = igdata.weights, bins = (xbins, ybins))

# ╔═╡ ceb02d26-9491-4fa7-b560-ff98a4733f45
mean(igdata.weights)

# ╔═╡ b0a81f9e-3900-410e-916e-f84b4e2494cb
plot(hgrad)

# ╔═╡ c8406ec4-1c96-45d3-ba89-a426aac7bac7
begin
hdump = deepcopy(h0)
umask = h0.weights .<= 0
hdump.weights[umask] .= 100
end

# ╔═╡ f629ac84-a3c8-4b7a-bc90-a8b6c6a999f2
plot(hdump)

# ╔═╡ 9d306ac1-085d-4270-b89a-9a33a4f97c43
plot(dhis[[0, 0]])

# ╔═╡ 100039de-14c0-4d89-ad5f-ea8258b752af
plot(dhis[[1, 0]])

# ╔═╡ ed041a49-9372-4e7a-80c3-2d883c9ed236
plot(dhis[[0, 0]])

# ╔═╡ c4991feb-a68b-40a0-b35b-71ee1255d04e
plot(ddelhis[[0, 0]])

# ╔═╡ 7c97dba0-a48e-45ae-aaad-c97fe3f58394
sum(ddelhis[[0, 0]].weights)

# ╔═╡ ef998fdc-eb6a-41d7-8a3a-f586467316df


# ╔═╡ 3aa4ab27-bd70-4dbb-bbd2-1b3cedb8e6a9
plot(hagrad)

# ╔═╡ 79ffe25f-a9d5-4097-8548-fa588a25e09f
for (i, move) in enumerate(moves)
	println(move, ", ", i, "\n ")
end

# ╔═╡ d5328a59-4ee8-48b3-83a6-648cb5bb394f
plot(dh)

# ╔═╡ 9cf31354-3425-4a84-938d-951a8efb5087
struct Frame
	steps::Tuple{Float64}
	bins::Array{Float64, 2}
	contents::Array{Float64, 2}
	mask::Array{Bool, 2}
end

# ╔═╡ edab54b4-941b-4750-8597-eac95485e0a1
md"""
## Compute deltas in each direction
"""

# ╔═╡ Cell order:
# ╠═5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╠═a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═1751666d-d26f-4d8e-b5d6-b36854dd7462
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╠═3922eba2-f322-4b06-b9e0-83bc723d7930
# ╠═3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╠═ceaff28e-ea37-4a92-a812-210fd2b91fad
# ╠═94bb1cff-3f57-4af6-b177-ff732cc78429
# ╠═40a8cad8-0a7b-41c5-9807-a5143c30a76b
# ╠═e4e0ac6b-ad55-4a96-81e7-c5ec8df53be6
# ╠═6faf5d14-f3f2-4549-bc1b-7d46771eb1da
# ╠═b05e2ea0-66c7-4c40-a998-92d7b8baa37f
# ╠═b7fb1cb1-941f-4bfa-ad96-7addfd1e60c0
# ╠═5cc3372f-f9df-4063-9e4c-dfa4e43f316c
# ╠═8f49762d-9e98-4fed-b55c-17d9425d0fff
# ╠═b9095b27-8f5d-40c5-9860-b7895a9955fb
# ╠═602dc3e2-bd48-4322-8687-42a03e1d0dd7
# ╠═b06f35cd-368a-4eec-b21d-68f70afe59fa
# ╠═df7b2cd3-5b2f-4b36-a200-aaedb38bd678
# ╠═62867032-c903-4b5b-b6db-51c02ddc6e6d
# ╠═33108e40-21a0-4ab7-bbc3-7827cb9c7a45
# ╠═c6b36bd0-97a1-4c02-8dcf-3d1f208299ba
# ╠═6f354604-3e53-4a6f-956b-eac04f01c949
# ╠═1d414c02-b275-4c27-a18f-99518f02442d
# ╠═50e73c5e-bc43-4cf6-b085-e1890e706e4a
# ╠═ce15bb7e-e9a1-4589-97b5-c209e4525f37
# ╠═98eabc22-3748-4b71-ae6b-44e98a6be83d
# ╠═211410f2-f516-4538-aa86-78aeac2d3b19
# ╠═5b0a5a5c-099d-45f4-ac0c-1c24bff169d0
# ╠═b5415946-1e85-432f-8a2e-c72eec552ee2
# ╠═6fc37709-d5ff-475f-877a-0b48bc54f02f
# ╠═f8f35d6b-fa66-4cd4-96ea-7e260fe2633e
# ╠═59c3d73d-2841-4340-abb0-4accfeb0bba9
# ╠═8c574ffb-dd95-4bb4-b0a3-e24d371b3c32
# ╠═1bf32c73-58db-46d2-b3d0-a6f2032ffafd
# ╠═69b2e4fb-7a76-4ddd-a172-0f2677b08648
# ╠═e204b639-cb6b-48d3-861b-83278bcd98b9
# ╠═465ccc2b-5d58-4a9a-9a56-24be0ceb971a
# ╠═d9c5642d-a93c-4094-b9b3-44b821ec20c8
# ╠═f524364b-1dd5-4445-b96d-b8f4c98504fd
# ╠═9da43e9a-d1cd-4ad5-9876-c41023f46fa0
# ╠═b5d78c84-1584-426e-9afa-c7d15297bedf
# ╠═d21157e9-7da5-449c-bb3a-de76084d0f62
# ╠═b049b708-934f-4531-89b8-66173c47d08c
# ╠═f8b2e785-50f4-4af2-b63c-85eeedba31ed
# ╠═6d5f5671-2c6f-460a-bb3b-4b1788f5efe4
# ╠═2339e4b0-e95b-4026-be77-e4306db471a5
# ╠═3930b27b-c1e8-4886-acbc-fc3e0099c5cb
# ╠═e9f7ff07-54c5-4f3a-a1d6-25d4d54499a2
# ╠═01cdd761-d68a-4563-914e-97e85b6438cb
# ╠═cfb9266e-cfc0-45ef-af6f-1ed5be7ce9f6
# ╠═8d50c759-262c-4e98-aa22-7e597a545e86
# ╠═f248c208-912a-40cb-9da5-c757af5301cc
# ╠═7f16716e-e452-4db8-8c6c-e600a0b50944
# ╠═2df93823-002f-46f7-9d64-b625a66fa688
# ╠═c22e56f9-1bab-499e-a80d-58afccc2ef8f
# ╠═a2f8b0b4-83ba-4c41-9a54-6f56d28febde
# ╠═eebf5cc7-58d0-44b4-9c2d-b65587bbff88
# ╠═fa808f07-96b6-4e1e-83e1-0694c969823c
# ╠═f9f106c6-14ff-46d0-b060-0a552c3a082c
# ╠═e12fb4a4-f24e-43db-9d12-deebc6f2b0ba
# ╠═c005b1d2-4ecf-4b44-bfcb-549ea746707d
# ╠═3b5730e5-145c-4fed-99b0-2e1da1982f68
# ╠═8f34afb7-f955-46d7-81ce-0eb282f23993
# ╠═ba0ba2a7-aa5b-4026-adde-3f767baafbe1
# ╠═b5c836bb-0155-41bd-88a7-70e62f2fcd3a
# ╠═5276e76f-5af5-49c3-80b7-b177f4eab65b
# ╠═eb9b5214-cfb9-4812-9e64-5b94ea419ca9
# ╠═a36ef6d9-b23d-4896-b930-11b026e762dc
# ╠═8c9d468a-ebf7-4873-b358-0aec1511cbe2
# ╠═0a7f02b0-e8b1-41fa-8bd9-741a0a1bcf11
# ╠═d174e2cd-d37f-48ac-8960-7dc2f815ddb5
# ╠═f91f8c7a-63e3-4518-bc8d-7bbd8c113c5b
# ╠═ceb02d26-9491-4fa7-b560-ff98a4733f45
# ╠═b0a81f9e-3900-410e-916e-f84b4e2494cb
# ╠═c8406ec4-1c96-45d3-ba89-a426aac7bac7
# ╠═f629ac84-a3c8-4b7a-bc90-a8b6c6a999f2
# ╠═9d306ac1-085d-4270-b89a-9a33a4f97c43
# ╠═100039de-14c0-4d89-ad5f-ea8258b752af
# ╠═ed041a49-9372-4e7a-80c3-2d883c9ed236
# ╠═c4991feb-a68b-40a0-b35b-71ee1255d04e
# ╠═7c97dba0-a48e-45ae-aaad-c97fe3f58394
# ╠═ef998fdc-eb6a-41d7-8a3a-f586467316df
# ╠═3aa4ab27-bd70-4dbb-bbd2-1b3cedb8e6a9
# ╠═79ffe25f-a9d5-4097-8548-fa588a25e09f
# ╠═d5328a59-4ee8-48b3-83a6-648cb5bb394f
# ╠═9cf31354-3425-4a84-938d-951a8efb5087
# ╠═edab54b4-941b-4750-8597-eac95485e0a1
