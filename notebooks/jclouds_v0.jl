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
import Graphs  as GG
import GraphPlot as GP 
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


# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""

## Generate Image

*data* is a namedtuple with the coordinates and energy
"""

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
img = fit(Histogram, (vec(xm), vec(ym)), weights(vec(zz)), (xbins, ybins))

# ╔═╡ 5cc3372f-f9df-4063-9e4c-dfa4e43f316c
plot(img)

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
data = histo_to_data(img)

# ╔═╡ 62867032-c903-4b5b-b6db-51c02ddc6e6d
histogram2d(data.x, data.y, weights = data.weights, bins = (xbins, ybins))

# ╔═╡ bc0a9ad4-62f4-4ef9-909a-4a2f68584b6f
md"""

## Gradient and Curvature functionality

**Internal functions**

* *_moves2d* define the movements in 2D (@TODO 3D)
* *_deltas* defines the deltas in each movement
* *_grad* computes the maximum gradient and its movement
* *_curvatures* computes the curvature for each movement
* *_minmax curvatures* compute the maximum and minimum curvature and their movements
* *_nodes* associate each cell to a node number
* *_neighbour node* associate each cell to the number of neighbour cells of other nodes, return the node number in any moves 

"""

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

# ╔═╡ cbf54f5a-1bdc-4473-b825-7d2f75211fc0
begin
	
function _deltas(coors, energy, edges, steps, m)
	
	xstep, ystep   = steps
	xx, yy         = coors
	xedges, yedges = edges
	
	his      = [fit(Histogram, (xx .- mx * xstep, yy .- my * ystep), weights(energy), (xedges, yedges)) for (mx, my) in m.moves]
	contents = deepcopy(his[m.imove0].weights)
	deltas   = [h.weights .- contents for h in his]
	
	dsteps   = [norm([xstep, ystep] .* move) for move in m.moves]
	dsteps[m.imove0] = 1.
	deltas  = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]
	
	return his, deltas
end

function _grad(deltas, m)
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


function _curvatures(deltas, m)

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
	
function _maxmin_curvatures(curves, m)
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

#--- Modes

function _node(i, j, igrad, m)
	index = CartesianIndex((i, j)...)
	imove = igrad[index]
	if (imove == m.imove0)
		return index
	else
		index = [i, j] + m.moves[imove]
		return _node(index[1], index[2], igrad, m)
	end
end

function _nodes(contents, igrad, m)

	cells = findall(x -> x.>0, contents)
	nodes = [_node(cell[1], cell[2], igrad, m) for cell in cells]
	unodes = unique(nodes)
	inodes = Dict()
	for (i, node) in enumerate(unodes)
		inodes[node] = i
	end
	
	xnodes = Int.(0 .* deepcopy(contents))
	for (k, cell) in enumerate(cells)
		xnodes[cell] = inodes[nodes[k]]
	end
	
	return xnodes
end

function _neighbour_node(coors, nodes, edges, steps, m)
	
	xstep, ystep   = steps
	xx, yy         = coors
	xedges, yedges = edges
	
	his      = [fit(Histogram, (vec(xx) .- mx * xstep, vec(yy) .- my * ystep), weights(vec(nodes)), (xedges, yedges)) for (mx, my) in m.moves]
	contents = deepcopy(his[m.imove0].weights)
	mask     = contents .> 0
	borders  = [(h.weights .>0) .* (h.weights .!= contents) for h in his] 
	nborders = reduce(.+, borders) .* mask
		
	return nborders, [h.weights for h in his]
end

#--- Public
	
function mesh(edges)
	nx, ny  = length(edges[1])-1, length(edges[2])-1
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]
	xm = ones(ny)'   .* centers[1]
	ym = centers[2]' .* ones(nx)
	return xm, ym
end
	
	
function quiver(idir, m, steps)
	xstep, ystep = steps
	uus = [m.dmoves[k] for k in vec(idir)]
	us = [xi * 0.8 * xstep for (xi, yi) in uus]
	vs = [yi * 0.8 * ystep for (xi, yi) in uus]	
	return us, vs
end
	
function kkclouds(coors, energy, steps)
	
	ndim  = length(coors)
	nsize = length(coors[1])

	# assert dimensions
	for i in 2:ndim
		@assert(length(coors[2]) == nsize)
	end
	@assert(length(energy) == nsize)
	@assert(length(steps)  == ndim)

	# define the extended edges 
	edges = [minimum(x)-1.5*step:step:maximum(x)+1.5*step 
		for (x, step) in zip(coors, steps)]

	# alias
	xx, yy         = coors
	xedges, yedges = edges
	
	# main histogram
	histo    = fit(Histogram, (xx, yy), weights(energy), (xedges, yedges))
	contents = deepcopy(histo.weights)

	# deltas
	m = _moves2d()
	his, deltas = _deltas(coors, energy, edges, steps, m)

	# gradient
	grad, igrad = _grad(deltas, m)

	# curvatures
	curves = _curvatures(deltas, m)
	lap    = curves[m.moves[m.imove0]]
	
	# maximum and monimum curvatures
	curmax, icurmax, curmin, icurmin = _maxmin_curvatures(curves, m)

	xs, ys = mesh(edges)
	
	# nodes
	xnodes = _nodes(contents, igrad, m)
	xborders, xneigh = _neighbour_node((xs, ys), xnodes, edges, steps, m)
	#nborders, neighbour_node = _neighbour_node(coors, xnodes, edges, steps, m)
	
	# output
	return (edges = (xedges, yedges),
		    x = xs, y = ys, contents = contents,
			grad = grad, igrad = igrad,
			lap = lap, curves = curves,
			curmax = curmax, icurmax = icurmax,
		    curmin = curmin, icurmin = icurmin,
			nodes  = xnodes,
			nborders = xborders, neighbour_node = xneigh)

end
	
end #begin

# ╔═╡ 441a1b38-7ff8-456c-8511-98f57828b26b
begin
	steps_   = [1., 1.]
	coors_   = (data.x, data.y)
	weights_ = data.weights
	mm       = _moves2d()
	cl       = kkclouds(coors_, weights_, steps_)
end;

# ╔═╡ 7e1adeb7-2714-4205-b895-62cbe0477008
keys(cl)

# ╔═╡ 7856646f-20f1-47b1-986b-0702e2f53305
md"""
### Control Histograms
"""

# ╔═╡ bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
cl.edges

# ╔═╡ b774e74d-0843-4f5c-98b9-9103bd0a5657
histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.contents), 
		bins = cl.edges, alpha = 0.8, title = "energy")
	

# ╔═╡ a684c917-3750-4d80-a86e-7cc3fd7b9d02
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.grad), 
		bins = cl.edges, alpha = 0.8, title = "gradient")
	quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.igrad, mm, steps_))
end

# ╔═╡ d8d579a5-aa12-4128-b2ff-252d165cd2a6
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.igrad), 
		bins = cl.edges, alpha = 0.8, title = "gradient direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.igrad, mm, steps_))
end

# ╔═╡ 1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.lap), 
		bins = cl.edges, alpha = 1., title = "laplacian")

# ╔═╡ fbdb0b59-f26f-42e5-a83a-a24962b09876
begin
mgrad = -2.0/sqrt(2.0)-1.0
md"""

Maximun grad $(mgrad)
"""
end

# ╔═╡ 55c5b8d3-9ac5-4ed1-8ec5-95cf74f81e6e
begin
chis2 = [histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.curves[m]), 
		bins = cl.edges, title = m) for m in mm.moves]
plot(chis2..., layout = (3, 3))
end

# ╔═╡ c7dbc705-c905-4fc8-b1a5-616345d029b8
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.curmax), 
		bins = cl.edges, alpha = 0.5, title = "curmax")
	quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmax, mm, steps_))
end

# ╔═╡ 62d11578-0d90-472d-8898-83e21c53d621
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.icurmax), 
		bins = cl.edges, alpha = 0.5, title = "curmax direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmax, mm, steps_))
end

# ╔═╡ 1605c6b4-f674-4209-ae46-f7ac4813693d
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.curmin), 
		bins = cl.edges, alpha = 0.5, title = "curmin direction")
	quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmin, mm, steps_))
end

# ╔═╡ f4155f75-af7f-4b79-b846-3bdc601d8767
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.icurmin), 
		bins = cl.edges, alpha = 0.5, title = "curmin direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmin, mm, steps_))
end

# ╔═╡ 654bce47-b51d-4c9a-81b3-b295f4bb055a
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.nodes), 
		bins = cl.edges, alpha = 1., title = "nodes")
end

# ╔═╡ 503b611f-d2f9-4b94-92c5-0b005505d5bf
begin
	histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.nborders), 
		bins = cl.edges, alpha = 1., title = "number of borders")
end

# ╔═╡ 27f27825-3f01-4c63-a119-4267ef69b11c
begin
nhis2 = [histogram2d(vec(cl.x), vec(cl.y), weights = vec(cl.neighbour_node[i]), 
		bins = cl.edges, title = m) for (i, m) in enumerate(mm.moves)]
plot(nhis2..., layout = (3, 3))
end

# ╔═╡ 713862ce-7003-4462-9e7b-d5611c3c96e2
md"""
## Graph

"""

# ╔═╡ e5d43a95-fcf1-4385-9e82-df1d2c70582c
begin
function _edges(nodes, neighs)
	iinodes = sort(unique(vec(nodes[nodes .>0])))
	dus     = Dict()
	#edges   = zeros(Bool, nnodes, nnodes)
	for iii in iinodes
		imask = cl.nodes .== iii
		us = []
		for neigh in neighs
			kmask = imask .* (neigh .> 0) .* (neigh .!= iii)
			ius   = unique(vec(neigh[kmask]))
			for k in ius
				if !(k in us)
					append!(us, k)
				end
			end
		end
		dus[iii] = us
		#for ui in us
		#	edges[iii, ui] = true
		#end
	end
	return dus#, edges
end

function cloud_graph(nodes, neighs)
	dus = _edges(nodes, neighs)

	nnodes = Base.size(unique(vec(cl.nodes)))[1] -1
	g = GG.Graph(nnodes)
	for inode in keys(dus)
		for knode in dus[inode]
			GG.add_edge!(g, inode, knode)
		end
	end
	return g
end
end # begin

# ╔═╡ ebd477fb-fef6-4ad9-8a25-c452172efa69
begin
gg = cloud_graph(cl.nodes, cl.neighbour_node)
GP.gplot(gg, nodelabel=1:GG.nv(gg), edgelabel=1:GG.ne(gg))
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
# ╠═7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╠═40a8cad8-0a7b-41c5-9807-a5143c30a76b
# ╠═e4e0ac6b-ad55-4a96-81e7-c5ec8df53be6
# ╠═6faf5d14-f3f2-4549-bc1b-7d46771eb1da
# ╠═b7fb1cb1-941f-4bfa-ad96-7addfd1e60c0
# ╠═5cc3372f-f9df-4063-9e4c-dfa4e43f316c
# ╠═8f49762d-9e98-4fed-b55c-17d9425d0fff
# ╠═602dc3e2-bd48-4322-8687-42a03e1d0dd7
# ╠═b06f35cd-368a-4eec-b21d-68f70afe59fa
# ╠═62867032-c903-4b5b-b6db-51c02ddc6e6d
# ╠═bc0a9ad4-62f4-4ef9-909a-4a2f68584b6f
# ╠═6890e500-b267-4005-859f-405ac0ef895a
# ╠═26bd9352-81d8-49de-a933-d62f22c462fa
# ╠═cbf54f5a-1bdc-4473-b825-7d2f75211fc0
# ╠═441a1b38-7ff8-456c-8511-98f57828b26b
# ╠═7e1adeb7-2714-4205-b895-62cbe0477008
# ╠═7856646f-20f1-47b1-986b-0702e2f53305
# ╠═bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
# ╠═b774e74d-0843-4f5c-98b9-9103bd0a5657
# ╠═a684c917-3750-4d80-a86e-7cc3fd7b9d02
# ╠═d8d579a5-aa12-4128-b2ff-252d165cd2a6
# ╠═1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
# ╠═fbdb0b59-f26f-42e5-a83a-a24962b09876
# ╠═55c5b8d3-9ac5-4ed1-8ec5-95cf74f81e6e
# ╠═c7dbc705-c905-4fc8-b1a5-616345d029b8
# ╠═62d11578-0d90-472d-8898-83e21c53d621
# ╠═1605c6b4-f674-4209-ae46-f7ac4813693d
# ╠═f4155f75-af7f-4b79-b846-3bdc601d8767
# ╠═654bce47-b51d-4c9a-81b3-b295f4bb055a
# ╠═503b611f-d2f9-4b94-92c5-0b005505d5bf
# ╠═27f27825-3f01-4c63-a119-4267ef69b11c
# ╠═713862ce-7003-4462-9e7b-d5611c3c96e2
# ╠═e5d43a95-fcf1-4385-9e82-df1d2c70582c
# ╠═ebd477fb-fef6-4ad9-8a25-c452172efa69
