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

# ╔═╡ b7fb1cb1-941f-4bfa-ad96-7addfd1e60c0
hh = fit(Histogram, (vec(xm), vec(ym)), weights(vec(zz)), (xbins, ybins))

# ╔═╡ 5cc3372f-f9df-4063-9e4c-dfa4e43f316c
plot(hh)

# ╔═╡ 8f49762d-9e98-4fed-b55c-17d9425d0fff
histogram2d(vec(xm), vec(ym), weights = vec(zz), nbins = (xbins, ybins))

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

# ╔═╡ 62867032-c903-4b5b-b6db-51c02ddc6e6d
histogram2d(data.x, data.y, weights = data.weights, bins = (xbins, ybins))

# ╔═╡ 6890e500-b267-4005-859f-405ac0ef895a
begin
	
struct Moves2D
	moves
	imove0
	dmoves
	omoves
end

function _moves2d()
	moves  = [[i, j] for i in -1:1:1 for j in -1:1:1]
	imove0 = [i for (i, move) in enumerate(moves) if move == [0, 0]][1]
	dmoves = Dict(1:9 .=> moves)
	omoves = Dict()
	imoves = deepcopy(moves)
	imoves = filter(move -> move != imoves[imove0], moves)
	for i in 1:9
		vals = 
		omoves[moves[i]] = [imove for imove in imoves 
			if dot(imove, moves[i]) == 0]	
	end
	str_moves = Moves2D(moves, imove0, dmoves, omoves)
	#str_moves.moves  = moves
	#str_moves.imove0 = imove0
	#str_moves.dmoces = dmoces
	#str_moves.omoces = omoves
end
	
end # begin

# ╔═╡ 26bd9352-81d8-49de-a933-d62f22c462fa
mov2d = _moves2d()

# ╔═╡ bcdee51d-3bfa-469b-b703-66a31762aebc
begin

function _jc_curvatures(deltas, m)

	ddeltas = Dict()
	for i in 1:9
		ddeltas[m.moves[i]] = deltas[i]
	end
	
	curvs = Dict()
	for imove in m.moves
		curvs[imove] = reduce(.+, [ddeltas[move] for move in m.omoves[imove]])
	end

	#dels = [reduce(.+, [ddeltas[i] for move in m.omoves[imove]]) for imove in m.moves] 

	return curvs
end
	
function _jc_maxmin_curvatures(curves, m)
	nx, ny = Base.size(curves[m.moves[m.imove0]])
	curmin  =  1e6 .* ones(nx, ny)
	icurmin = m.imove0 .* ones(Int, nx, ny) 
	curmax  = -1e6 .* ones(nx, ny)
	icurmax = m.imove0 .* ones(Int, nx, ny) 
	for (i, move) in enumerate(m.moves)
		dd = curves[move]
		mask1 = dd .> curmax
		curmax[mask1] .= dd[mask1]
		icurmax[mask1] .= i
		mask2 = dd .< curmin
		curmin[mask2] .= dd[mask2]
		icurmin[mask2] .= i 
	end
	return curmax, icurmax, curmin, icurmin
end

end # begin

# ╔═╡ cbf54f5a-1bdc-4473-b825-7d2f75211fc0
begin

struct JClouds
	steps
	edges
	centers
	weights
	deltas
	curves
	grad
	igrad
	curmax
	icurmax
	curmin
	icurmin
end

function _jc_deltas(coors, energy, edges, steps, m)
	
	xstep, ystep   = steps
	xx, yy         = coors
	xedges, yedges = edges
	
	his    = [fit(Histogram, (xx .- mx * xstep, yy .- my * ystep), weights(energy), (xedges, yedges)) for (mx, my) in m.moves]
	contents = deepcopy(his[m.imove0].weights)
	deltas   = [h.weights .- contents for h in his]
	
	dsteps   = [norm([xstep, ystep] .* move) for move in m.moves]
	dsteps[m.imove0] = 1.
	deltas  = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]
	
	return his, deltas
end

function _jc_grad(deltas, m)
	nx, ny = Base.size(deltas[m.imove0])
	d0     = deltas[m.imove0]
	grad   = deepcopy(d0)
	igrad  = m.imove0 .* ones(Int, nx, ny)
	for (i, di) in enumerate(deltas)
		imask         = di .> grad
		grad[imask]  .= di[imask]
		igrad[imask] .= i
	end
	return grad, igrad
end


	
function _jc_mesh(edges)
	nx, ny  = length(edges[1])-1, length(edges[2])-1
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]
	xm = ones(ny)'   .* centers[1]
	ym = centers[2]' .* ones(nx)
	return xm, ym
end
	
	
function _jc_quiver(idir, m, steps)
	xstep, ystep = steps
	uus = [m.dmoves[k] for k in vec(idir)]
	us = [xi * 0.8 * xstep for (xi, yi) in uus]
	vs = [yi * 0.8 * ystep for (xi, yi) in uus]	
	return us, vs
end
	
function _jclouds(coors, energy, steps)
	ndim = length(coors)
	nsize = length(coors[1])

	# assert dimensions
	for i in 2:ndim
		@assert(length(coors[2]) == nsize)
	end
	@assert(length(energy)  == nsize)
	@assert(length(steps) == ndim)

	# create main histogram
	edges = [minimum(x)-1.5*step:step:maximum(x)+1.5*step 
		for (x, step) in zip(coors, steps)]
	#his    = [fit(Histogram, (data.x .- mx * xstep, data.y .- my * ystep), #weights(data.weights), (xbins, ybins)) for (mx, my) in moves]

	# alias
	xx, yy         = coors
	xedges, yedges = edges
	
	# main histogram
	histo = fit(Histogram, (xx, yy), weights(energy), (xedges, yedges))
	contents = deepcopy(histo.weights)

	# deltas
	m = _moves2d()
	his, deltas = _jc_deltas(coors, energy, edges, steps, m)

	# gradient
	grad, igrad = _jc_grad(deltas, m)

	# curvatures
	curves = _jc_curvatures(deltas, m)
	lap    = curves[m.moves[m.imove0]]
	
	# maximum and monimum curvatures
	curmax, icurmax, curmin, icurmin = _jc_maxmin_curvatures(curves, m)
	
	
	return (edges = edges, contents = contents,
			grad = grad, igrad = igrad,
			lap = lap,
			curmax = curmax, icurmax = icurmax,
		    curmin = curmin, icurmin = icurmin)
	
end
	
end #begin

# ╔═╡ 441a1b38-7ff8-456c-8511-98f57828b26b
begin
	steps_   = [1., 1.]
	coors_   = (data.x, data.y)
	weights_ = data.weights
	mm       = _moves2d()
	jclouds  = _jclouds(coors_, weights_, steps_)
end

# ╔═╡ bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
jclouds.edges

# ╔═╡ 21593ffc-e07f-4d67-a683-b7e87e42f086
begin
	xxx, yyy = _jc_mesh(jclouds.edges)
	xxedges  = jclouds.edges[1], jclouds.edges[2]
	histogram2d(vec(xxx), vec(yyy), weights = vec(jclouds.grad), 
		bins = xxedges, alpha = 0.5, title = "gradient")
	quiver!(vec(xxx), vec(yyy), quiver = _jc_quiver(jclouds.igrad, mm, steps_))
end

# ╔═╡ 1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
histogram2d(vec(xxx), vec(yyy), weights = vec(jclouds.curmin), 
		bins = xxedges, alpha = 0.5, title = "laplacian")

# ╔═╡ c7dbc705-c905-4fc8-b1a5-616345d029b8
begin
	histogram2d(vec(xxx), vec(yyy), weights = vec(jclouds.curmax), 
		bins = xxedges, alpha = 0.5, title = "curmax direction")
	quiver!(vec(xxx), vec(yyy), quiver = _jc_quiver(jclouds.icurmax, mm, steps_))
end

# ╔═╡ f83362d1-f485-4ea9-a75a-32424bf9dbbb
begin
	histogram2d(vec(xxx), vec(yyy), weights = vec(jclouds.curmin), 
		bins = xxedges, alpha = 0.5, title = "curmin direction")
	quiver!(vec(xxx), vec(yyy), quiver = _jc_quiver(jclouds.icurmin, mm, steps_))
end

# ╔═╡ 713862ce-7003-4462-9e7b-d5611c3c96e2
md"""
### Debugging
"""

# ╔═╡ e12fb4a4-f24e-43db-9d12-deebc6f2b0ba
moves[4]

# ╔═╡ 3b5730e5-145c-4fed-99b0-2e1da1982f68
moves[1], moves[2]

# ╔═╡ b5c836bb-0155-41bd-88a7-70e62f2fcd3a
moves[1]

# ╔═╡ 79ffe25f-a9d5-4097-8548-fa588a25e09f
for (i, move) in enumerate(moves)
	println(move, ", ", i, "\n ")
end

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
# ╠═b7fb1cb1-941f-4bfa-ad96-7addfd1e60c0
# ╠═5cc3372f-f9df-4063-9e4c-dfa4e43f316c
# ╠═8f49762d-9e98-4fed-b55c-17d9425d0fff
# ╠═602dc3e2-bd48-4322-8687-42a03e1d0dd7
# ╠═b06f35cd-368a-4eec-b21d-68f70afe59fa
# ╠═62867032-c903-4b5b-b6db-51c02ddc6e6d
# ╠═6890e500-b267-4005-859f-405ac0ef895a
# ╠═26bd9352-81d8-49de-a933-d62f22c462fa
# ╠═bcdee51d-3bfa-469b-b703-66a31762aebc
# ╠═cbf54f5a-1bdc-4473-b825-7d2f75211fc0
# ╠═441a1b38-7ff8-456c-8511-98f57828b26b
# ╠═bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
# ╠═21593ffc-e07f-4d67-a683-b7e87e42f086
# ╠═1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
# ╠═c7dbc705-c905-4fc8-b1a5-616345d029b8
# ╠═f83362d1-f485-4ea9-a75a-32424bf9dbbb
# ╠═713862ce-7003-4462-9e7b-d5611c3c96e2
# ╠═e12fb4a4-f24e-43db-9d12-deebc6f2b0ba
# ╠═3b5730e5-145c-4fed-99b0-2e1da1982f68
# ╠═b5c836bb-0155-41bd-88a7-70e62f2fcd3a
# ╠═79ffe25f-a9d5-4097-8548-fa588a25e09f
