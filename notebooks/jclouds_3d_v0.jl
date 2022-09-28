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

Santiago, August 2022

---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly()

# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""

## Generate Image

A small img 2D cataloge to test clouds 2D

"""

# ╔═╡ d8a02e0a-db35-4965-8322-8741c3ffbd49
begin
"""
Just a smeared line!
"""
function line(; ndim = 3, threshold = 0.)

	tstep = 0.1
	ts = 0:tstep:1.
	ax, bx, cx = 5., 0., 1.
	ay, by, cy = -5., -5., -1.
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

# ╔═╡ 2be6a8f8-74f2-49d3-95b1-6add8103c310
begin

#imgs = Dict()
#imgs[:line]  = line()
#imgs[:braid] = braid()

#bimg = @bind keyimg Select([:line, :braid])

md"""

Selecting line

"""
end

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin
bndim = @bind nndim Select([2, 3])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select label to plot $(bndim)

"""
end

# ╔═╡ 6c8bf138-8fec-4c69-b4dd-4284faddeed0
begin
img = line(ndim = nndim, threshold = 6.)
end;

# ╔═╡ 5a1832c1-33ff-45dc-8f47-212179dbe862
md"""

## Clouds code

  * Moves 
"""

# ╔═╡ 9c413a7e-3cbc-4676-b846-37db33f0c2f0
begin

struct Moves
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
		omoves[moves[i]] = [imove for imove in imoves
			if LA.dot(imove, moves[i]) == 0]
	end
	smoves = Moves(moves, imove0, dmoves, omoves)
	return smoves
end

function _moves3d()
	moves  = [[i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1]
	nmoves = length(moves)
	imove0 = [i for (i, move) in enumerate(moves) if move == [0, 0, 0]][1]
	dmoves = Dict(1:nmoves .=> moves)
	omoves = Dict()
	imoves = deepcopy(moves)
	imoves = filter(move -> move != imoves[imove0], moves)
	for i in 1:nmoves
		omoves[moves[i]] = [imove for imove in imoves
			if LA.dot(imove, moves[i]) == 0]
	end
	smoves = Moves(moves, imove0, dmoves, omoves)
	return smoves
end

end

# ╔═╡ 1367811b-355c-446d-9d58-d07a71dfdb23
begin

function _hcoors(ucoors, move)
	ndim = length(move)
	z    = ucoors .- move'
	zt   = Tuple(z[:, i] for i in 1:ndim)
	return zt
end

function _deltas(ucoors, energy, edges, steps, m)

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move), 
		SB.weights(energy), edges) for move in m.moves]
    contents = deepcopy(his[m.imove0].weights)
    deltas   = [h.weights .- contents for h in his]

    dsteps   = [LA.norm(steps .* move) for move in m.moves]
    dsteps[m.imove0] = 1.
    deltas  = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]

    return deltas
end

function _gradient(deltas, m)
    dims   = Base.size(deltas[m.imove0])
    d0     = deltas[m.imove0]
    grad   = deepcopy(d0)
    igrad  = m.imove0 .* ones(Int, dims...)
    for (i, di) in enumerate(deltas)
        imask         = di .> grad
        grad[imask]  .= di[imask]
        igrad[imask] .= i
    end
    return grad, igrad
end


function _curvatures(deltas, m)

    ddeltas = Dict()
    for i in 1:length(m.moves)
        ddeltas[m.moves[i]] = deltas[i]
    end
    curvs = Dict()
    for imove in m.moves
        curvs[imove] = reduce(.+, [ddeltas[move] for move in m.omoves[imove]])
    end
    return curvs
end

function _maxmin_curvatures(curves, m)
    nsize = Base.size(curves[m.moves[m.imove0]])
    curmin  =  1e6 .* ones(nsize...)
    icurmin = m.imove0 .* ones(Int, nsize...)
    curmax  = -1e6 .* ones(nsize...)
    icurmax = m.imove0 .* ones(Int, nsize...)
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

function _node(cell, igrad, m)
    imove = igrad[cell]
    if (imove == m.imove0)
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
	
    contents = deepcopy(his[m.imove0].weights)
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

# ╔═╡ 23c6a3f5-443f-4199-87c8-68b93e495107
begin
	
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
	moves  = ndim == 2 ? _moves2d() : _moves3d()
	move0  = moves.moves[moves.imove0]
	ucoors = reduce(hcat, coors)

    # main histogram
    histo    = SB.fit(SB.Histogram, _hcoors(ucoors, move0), SB.weights(energy), edges)
    contents = deepcopy(histo.weights)
	cells    = findall(x -> x .> threshold, contents)

    # deltas
    deltas = _deltas(ucoors, energy, edges, steps, moves)

    # gradient
    grad, igrad = _gradient(deltas, moves)

    # curvatures
    curves = _curvatures(deltas, moves)
    lap    = curves[move0]

    # maximum and monimum curvatures
    curmax, icurmax, curmin, icurmin = _maxmin_curvatures(curves, moves)

    #xs, ys = mesh(edges)

    # nodes
    xnodes = _nodes(igrad, cells, moves)
    xborders, xneigh = _neighbour_node(ucoors, xnodes, edges, steps, moves)

	#xlinks = _links(xnodes, xneigh)
	
    # output
    return (edges = edges,  coors = coors, contents = contents[cells],
			cells = cells,
            grad = grad[cells], igrad = igrad[cells],
            lap = lap[cells], curves = curves,
            curmax = curmax[cells], icurmax = icurmax[cells],
            curmin = curmin[cells], icurmin = icurmin[cells],
            nodes  = xnodes,
            nborders = xborders[cells], neighbour_node = xneigh)
#			links = xlinks)

end
end # begin

# ╔═╡ d2c45f38-8cbe-474e-b66a-e45cfaf72f50
md"""

# ╔═╡ f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
begin
steps_ = [edge[2]-edge[1] for edge in img.edges]
xcl   = clouds(img.coors, img.contents, steps_)
end;

# ╔═╡ 13ac9fdf-46d0-4940-80e3-8619f0609108
md"""

## Plots
"""

# ╔═╡ a689debb-8763-45c4-a03d-94c8e970b243
begin

blabel = @bind label Select([:contents, :grad, :igrad, :lap, :curmax, :icurmax, :curmin, :icurmin, :nodes, :nborders])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select label to plot $(blabel)

"""
end

# ╔═╡ dfa64554-5fb1-4d63-80d3-19aee7a476b8
begin
function cplot(cl, label, title)
	ndim = length(cl.coors)
	vals = getfield(cl, label)
	p1 = ndim == 2 ? histogram2d(cl.coors..., weights = vals, nbins = cl.edges) : p1 = scatter(cl.coors..., zcolor = vals, alpha = 0.1)
	p2 = histogram(vals, nbins = 100)
	plot(p1, p2, title = title)
end
end

# ╔═╡ fe2513f9-461c-4066-9e30-6c55540ccd4b
cplot(img, :contents, :contents)

# ╔═╡ 1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
cplot(xcl, label, label)

# ╔═╡ 16b988b0-887f-4672-b347-9c374fcc3fae
begin

function cloud_graph(nodes, neighs)
	dus = _links(nodes, neighs)

	nnodes = length(unique(nodes))
	g = GG.Graph(nnodes)
	for inode in keys(dus)
		for knode in dus[inode]
			GG.add_edge!(g, inode, knode)
		end
	end
	return g
end
end # begin

# ╔═╡ 1c402508-afd3-46a1-8dbc-a23fd9bd63e1
begin
gg = cloud_graph(xcl.nodes, xcl.neighbour_node)
GP.gplot(gg, nodelabel=1:GG.nv(gg), edgelabel=1:GG.ne(gg))
end

# ╔═╡ 068e9533-1e4a-40be-83db-617a17935b0c
"""
linear scale a vector of values between a minimum and a maximum

Parameters:
	var : Vector{Real}
	emin: Real
	emax: Real

Return:
	vvar: Vector{Real}

"""
function vscale(var, emin = 1, emax = 4)
	vvar = (var .- minimum(var)) ./(maximum(var) - minimum(var)) .*(emax-emin) .+ emin
	return vvar
end

# ╔═╡ 654fb60f-2349-4b22-934e-dfb43080f5ec
begin
function _mesh3d(x, y, z)
	xv = getindex.(Iterators.product(x, y, z), 1)  # first.(Iterators.product(x, y, z), 1) is also ok
	yv = getindex.(Iterators.product(x, y, z), 2)
	zv = getindex.(Iterators.product(x, y, z), 3)
	return xv, yv, zv
end
end #begin

# ╔═╡ Cell order:
# ╠═5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╠═a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╠═3922eba2-f322-4b06-b9e0-83bc723d7930
# ╠═3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╠═7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╠═d8a02e0a-db35-4965-8322-8741c3ffbd49
# ╠═2be6a8f8-74f2-49d3-95b1-6add8103c310
# ╠═e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╠═6c8bf138-8fec-4c69-b4dd-4284faddeed0
# ╠═fe2513f9-461c-4066-9e30-6c55540ccd4b
# ╠═5a1832c1-33ff-45dc-8f47-212179dbe862
# ╠═9c413a7e-3cbc-4676-b846-37db33f0c2f0
# ╠═1367811b-355c-446d-9d58-d07a71dfdb23
# ╠═23c6a3f5-443f-4199-87c8-68b93e495107
# ╠═d2c45f38-8cbe-474e-b66a-e45cfaf72f50
# ╠═f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
# ╠═13ac9fdf-46d0-4940-80e3-8619f0609108
# ╠═a689debb-8763-45c4-a03d-94c8e970b243
# ╟─dfa64554-5fb1-4d63-80d3-19aee7a476b8
# ╠═1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
# ╠═16b988b0-887f-4672-b347-9c374fcc3fae
# ╠═1c402508-afd3-46a1-8dbc-a23fd9bd63e1
# ╠═068e9533-1e4a-40be-83db-617a17935b0c
# ╠═654fb60f-2349-4b22-934e-dfb43080f5ec
