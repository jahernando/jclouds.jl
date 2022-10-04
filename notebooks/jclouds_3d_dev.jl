### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5dcb2929-115e-459c-b98d-43ae7bcabd3a
using Pkg; Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/jclouds")

# ╔═╡ a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
begin
using Markdown
#using HDF5
#using DataFrames
import StatsBase as SB
using Plots
import LinearAlgebra as LA
#using Statistics
using PlutoUI
using Random
using Base
import Images.ImageFiltering as IF
#using TestImages
#using Distributions
import Graphs  as GG
import GraphPlot as GP
end

# ╔═╡ a57cdb41-c388-4976-bec8-ec0650fb139c
import jclouds as jc

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


Dev NB for clouds 3D in Julia


J.A. Hernado,

Santiago, September 2022

---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly();

# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""

## Generate Image

Produces a smeared curve in 2D or ·D

"""

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin
bndim = @bind nndim Select([2, 3])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select dimensions of the line $(bndim)

"""
end

# ╔═╡ 1a8e9aa9-a47d-40fd-84c6-cfa49f9b1cc4
begin

md"""

**Image**

"""
end

# ╔═╡ 5a1832c1-33ff-45dc-8f47-212179dbe862
md"""

## Clouds 
"""

# ╔═╡ 13ac9fdf-46d0-4940-80e3-8619f0609108
md"""

## Plots
"""

# ╔═╡ a689debb-8763-45c4-a03d-94c8e970b243
begin

blabel = @bind label Select([:contents, :grad, :igrad, :lap, :curmax, :icurmax, :curmin, :icurmin, :nodes, :nborders])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select label to plot :  $(blabel)

"""
end

# ╔═╡ 7b7981ca-1540-48a1-88e1-4f27e7787b70
md"""
## Graph
"""

# ╔═╡ 4751c9f4-e6d2-4b4a-b34f-7dbb258d443a
md"""

## DEV Code


"""

# ╔═╡ b3d9c017-f1ae-4ec8-8bad-8adf1774a7c7
md"""
## Extra code
"""

# ╔═╡ 501d3779-c536-4ae4-b046-5cbe95fab4e5
begin

#----------
# Moves
#----------

"""
Function that returns a NamedTuple with the one step movement vectors
in 2D and 3D.

It also returns the null move, the symetric pair of movees,
and the ortogonal move associated to any symmetric movement.
The null (*i0*), symmetric (*isym*), and ortogonal movs (*iorto*)
are indexed respect the main vector of moves (*moves*)

i.e 2D: null [0., 0.], symmetric pair ([1, 0], [-1, 0]),
ortogonal directions of the previous pair ([0, 1], [0, -1])

"""
function moves(ndim)
	moves = ndim == 2 ? [[i, j] for i in -1:1:1 for j in -1:1:1] : [[i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1]
	move0 = ndim == 2 ? [0, 0] : [0, 0, 0]
	kmove0 = [i for (i, move) in enumerate(moves) if (move == move0)][1]
	smoves = [(i, j) for (i, movei) in enumerate(moves) for (j, movej) in enumerate(moves) if (movei == -1 .* movej ) &  (i > j)]
	omoves = Dict()
	for ii in smoves
		movei = moves[ii[1]]
		omoves[ii] = [j for (j, movej) in enumerate(moves) if ((sum(movej .* movei) == 0.0) & (sum(movej .* movej) != 0.0))]
	end
	return (moves = moves, i0 = kmove0, isym = smoves, iortho = omoves)
end

#------------
# Clouds
#------------


#------------
# Clouds
#------------

"""

Create Clouds

From a (x, y) or (x, y, z) tuples of points with an energy.
The space is binned in steps and a histogram is created with the energy.
Only bins with contents > threshold are considered valid cells.

For each cell several information is provided:

content, gradient, laplacian, minimum and maximum curvature, node number,
and number of neighbour cells belonging to another node.

A dictionary, nodes_edges, with the edes between the nodes is also provided.

"""
function clouds(coors, energy, steps, threshold = 0.)

    ndim  = length(coors)
    nsize = length(coors[1])

    # assert dimensions
    for i in 2:ndim
        @assert(length(coors[i]) == nsize)
    end
    @assert(length(energy) == nsize)
    @assert(length(steps)  == ndim)

    # define the extended edges
    edges = Tuple(minimum(x) - 1.5*step : step : maximum(x) + 1.5*step for (x, step) in zip(coors, steps))

    # alias
	mm  = moves(ndim)
	m0  = mm.moves[mm.i0]
	ucoors = reduce(hcat, coors)

    # main histogram
    histo    = SB.fit(SB.Histogram, _hcoors(ucoors, m0), SB.weights(energy), edges)
    contents = deepcopy(histo.weights)
	cells    = findall(x -> x .> threshold, contents)

    # deltas
    deltas = _deltas(ucoors, energy, edges, steps, mm, threshold)

    # gradient
    grad, igrad = _gradient(deltas, mm)

    # curvatures
    curves = _curvatures(deltas, mm)
    lap    = reduce(.+, deltas)

    # maximum and monimum curvatures
    curmax, icurmax, curmin, icurmin = _maxmin_curvatures(Base.size(lap), curves, mm)

    #xs, ys = mesh(edges)

    # nodes
    xnodes = _nodes(igrad, cells, mm)
    xborders, xneigh = _neighbour_node(ucoors, xnodes, edges, steps, mm)

	xlinks = _links(xnodes, xneigh)

    # output
    return (edges = edges,  coors = coors, contents = contents[cells],
			cells = cells,
            grad = grad[cells], igrad = igrad[cells],
            lap = lap[cells], #curves = curves,
            curmax = curmax[cells], icurmax = icurmax[cells],
            curmin = curmin[cells], icurmin = icurmin[cells],
            nodes  = xnodes,
            nborders = xborders[cells], nodes_edges = xlinks)

end


#-----------------------------
#  Clouds internal functions
#-----------------------------

function _hcoors(ucoors, move)
	ndim = length(move)
	z    = ucoors .- move'
	zt   = Tuple(z[:, i] for i in 1:ndim)
	return zt
end

function _deltas(ucoors, energy, edges, steps, m, threshold)

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
		SB.weights(energy), edges) for move in m.moves]
    contents = deepcopy(his[m.i0].weights)
	asize = Base.size(contents)
	mask = contents .<= threshold
	deltas   = [h.weights .- contents for h in his]
	for (delta, h) in zip(deltas, his)
		delta[h.weights .<= threshold] .= 0.0
		delta[mask] .= 0.0
	end
	
    dsteps   = [LA.norm(steps .* move) for move in m.moves]
    dsteps[m.i0] = 1.
    deltas  = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]

    return deltas
end

function _gradient(deltas, m)
    dims   = Base.size(deltas[m.i0])
    d0     = deltas[m.i0]
    grad   = deepcopy(d0)
    igrad  = m.i0 .* ones(Int, dims...)
    for (i, di) in enumerate(deltas)
        imask         = di .> grad
        grad[imask]  .= di[imask]
        igrad[imask] .= i
    end
    return grad, igrad
end


function _curvatures(deltas, m)

    curvs = Dict()
    for smove in m.isym
       curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in m.iortho[smove]])
		#curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in smove])
    end
    return curvs
end

function _maxmin_curvatures(nsize, curves, m)
    curmin  =  1e6 .* ones(nsize...)
    icurmin = m.i0 .* ones(Int, nsize...)
    curmax  = -1e6 .* ones(nsize...)
    icurmax = m.i0 .* ones(Int, nsize...)
    for smove in m.isym
		imove = smove[1]
        dd = curves[imove]
        mask1 = dd .> curmax
        curmax[mask1] .= dd[mask1]
        icurmax[mask1] .= imove
        mask2 = dd .< curmin
        curmin[mask2] .= dd[mask2]
        icurmin[mask2] .= imove
	end
    return curmax, icurmax, curmin, icurmin
end

function _node(cell, igrad, m)
    imove = igrad[cell]
    if (imove == m.i0)
        return cell
    else
        #cindex_ = tuple(cindex) + m.moves[imove]
		nextcell = Tuple(Tuple(cell) .+ m.moves[imove])
        return _node(CartesianIndex(nextcell), igrad, m)
    end
end

function _nodes(igrad, cells, m)

	cnodes  = [_node(cell, igrad, m) for cell in cells]
	ucnodes = unique(cnodes)
	dicnodes = Dict()
	for (i, cnode) in enumerate(ucnodes)
    	dicnodes[Tuple(cnode)] = i
	end
	nodes = [dicnodes[Tuple(cnode)] for cnode in cnodes]
	return nodes
end

function _neighbour_node(ucoors, nodes, edges, steps, m)

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
		SB.weights(nodes), edges) for move in m.moves]

    contents = deepcopy(his[m.i0].weights)
    mask     = contents .> 0
    borders  = [(h.weights .>0) .* (h.weights .!= contents) for h in his]
    nborders = reduce(.+, borders) .* mask

    return nborders, [h.weights for h in his]
end


function _links(nodes, neighs)

	imove0 = length(neighs) > 9 ? 14 : 5

	dus     = Dict()
	for inode in nodes
		imask = neighs[imove0] .== inode
		us = []
		for neigh in neighs
			kmask = imask .* (neigh .> 0) .* (neigh .!= inode)
			ius   = unique(vec(neigh[kmask]))
			for k in ius
				if !(k in us)
					append!(us, k)
				end
			end
		end
		dus[inode] = us
	end
	return dus
end

end

#

# ╔═╡ d8a02e0a-db35-4965-8322-8741c3ffbd49
begin
"""
Just a smeared line!
"""
function line(;ndim = 3, threshold = 0.)

	tstep = 0.1
	ts = 0:tstep:1.
	ax, bx, cx = 5., 0., 0.
	ay, by, cy = -5., -5., 0.
	az, bz, cz = 1., -1., 0.
	xx = cx .* ts .* ts + ax .* ts .+ bx
	yy = cy .* ts .* ts + ay .* ts .+ by
	zz = cz .* ts .* ts + az .* ts .+ bz 

	zsig  = 5.
	sigma = 2 * tstep #* ma(ax, ay, az)
	xxbins = minimum(xx) - zsig .* sigma : sigma : maximum(xx) + zsig .*sigma
	yybins = minimum(yy) - zsig .* sigma : sigma : maximum(yy) + zsig .*sigma
	zzbins = minimum(zz) - zsig .* sigma : sigma : maximum(zz) + zsig .*sigma 

	heigth  = 1000
	xxcontents = heigth * ones(Base.size(xx))

	coors = ndim == 2 ? (xx, yy) : (xx, yy, zz)
	edges = ndim == 2 ? (xxbins, yybins) : (xxbins, yybins, zzbins)
	hh  = SB.fit(SB.Histogram, coors, SB.weights(xxcontents), edges)

	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]

	factor = 10.
	sigma_kernel = (sigma * factor) .* ones(ndim)
	weights_ = IF.imfilter(hh.weights, IF.Kernel.gaussian(sigma_kernel))

	cells    = findall(x -> x .> threshold, weights_)
	coors    = [[centers[i][cell[i]] for cell in cells] for i in 1:1:ndim]
	contents = [hh.weights[cell] for cell in cells]

	contents = [weights_[cell] for cell in cells]

	return (coors = coors, cells = cells, contents = contents, edges = edges)

end

end # begin

# ╔═╡ 6c8bf138-8fec-4c69-b4dd-4284faddeed0
begin
img = line(ndim = nndim, threshold = 6.)
end;

# ╔═╡ f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
begin
steps = [edge[2]-edge[1] for edge in img.edges]
xcl   = jc.clouds(img.coors, img.contents, steps)
end;

# ╔═╡ 4e43c8e3-89e2-44ca-a6ed-48a364d90486
begin
md"""

steps of the voxels: $(steps[1])

"""
end

# ╔═╡ 0e47793c-3623-4033-89f8-c5b2e90e9c5b
begin
vals = getfield(xcl, label)
minv, maxv = minimum(vals), maximum(vals)
end;

# ╔═╡ 3c20ca80-2ad6-42b2-9612-3f6ec25bc21a
begin
brange0 = @bind v0 Slider(minv:maxv, default = minv)
brange1 = @bind v1 Slider(minv:maxv, default = maxv)
md"""

Selec range for variable $(label):

minimum $(brange0)
maximum  $(brange1)

"""
end

# ╔═╡ c0e4ca98-e666-46b6-a693-0717ea39fad0
md"""

Selected range : [ $(v0), $(v1) ]

"""

# ╔═╡ dfa64554-5fb1-4d63-80d3-19aee7a476b8
begin
function cplot(cl, label, title, vrange)
	ndim = length(cl.coors)
	vals = getfield(cl, label)
	mask = (vals .>= vrange[1]) .* (vals .<= vrange[2])
	coors = [c[mask] for c in cl.coors]
	vvals = vals[mask]
	theme(:dark)
	p1 = ndim == 2 ? histogram2d(coors..., weights = vvals, nbins = cl.edges) : p1 = scatter(coors..., zcolor = vvals, alpha = 0.1)
	p2 = histogram(vvals, nbins = 100)
	plot(p1, p2, title = title)
end
end

# ╔═╡ d26c89ae-1629-4e98-8bde-3e8abe8bfd8d
cplot(img, :contents, :contents, [minimum(img.contents), maximum(img.contents)])

# ╔═╡ 1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
cplot(xcl, label, label, [v0, v1])

# ╔═╡ 16b988b0-887f-4672-b347-9c374fcc3fae
begin

function cloud_graph(nodes, nodes_edges)
	
	nnodes = length(unique(nodes))
	g = GG.Graph(nnodes)
	for inode in keys(nodes_edges)
		for knode in nodes_edges[inode]
			GG.add_edge!(g, inode, knode)
		end
	end
	return g
end
	
end # begin

# ╔═╡ 1c402508-afd3-46a1-8dbc-a23fd9bd63e1
begin
gg = cloud_graph(xcl.nodes, xcl.nodes_edges)
GP.gplot(gg, nodelabel=1:GG.nv(gg), edgelabel=1:GG.ne(gg))
end

# ╔═╡ Cell order:
# ╟─5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╟─a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─3922eba2-f322-4b06-b9e0-83bc723d7930
# ╟─3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╟─7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╟─e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╠═6c8bf138-8fec-4c69-b4dd-4284faddeed0
# ╟─1a8e9aa9-a47d-40fd-84c6-cfa49f9b1cc4
# ╠═d26c89ae-1629-4e98-8bde-3e8abe8bfd8d
# ╟─5a1832c1-33ff-45dc-8f47-212179dbe862
# ╠═f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
# ╟─4e43c8e3-89e2-44ca-a6ed-48a364d90486
# ╟─13ac9fdf-46d0-4940-80e3-8619f0609108
# ╟─a689debb-8763-45c4-a03d-94c8e970b243
# ╟─0e47793c-3623-4033-89f8-c5b2e90e9c5b
# ╟─3c20ca80-2ad6-42b2-9612-3f6ec25bc21a
# ╟─c0e4ca98-e666-46b6-a693-0717ea39fad0
# ╠═1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
# ╟─7b7981ca-1540-48a1-88e1-4f27e7787b70
# ╟─1c402508-afd3-46a1-8dbc-a23fd9bd63e1
# ╠═4751c9f4-e6d2-4b4a-b34f-7dbb258d443a
# ╟─b3d9c017-f1ae-4ec8-8bad-8adf1774a7c7
# ╠═501d3779-c536-4ae4-b046-5cbe95fab4e5
# ╠═d8a02e0a-db35-4965-8322-8741c3ffbd49
# ╠═dfa64554-5fb1-4d63-80d3-19aee7a476b8
# ╠═16b988b0-887f-4672-b347-9c374fcc3fae
