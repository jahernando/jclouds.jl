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

	sigma = 0.1
	ts = 0:0.05:1.
	ax, bx, cx = 1., 0., 0.
	ay, by, cy = -1., -5., 0.
	az, bz, cz = 0.5, -1., 0.
	xx = cx .* ts .* ts + ax .* ts .+ bx
	yy = cy .* ts .* ts + ay .* ts .+ by
	zz = cz .* ts .* ts + az .* ts .+ bz


	zsig   = 6
	xxbins = minimum(xx) - zsig .* sigma : sigma : maximum(xx) + zsig .*sigma
	yybins = minimum(yy) - zsig .* sigma : sigma : maximum(yy) + zsig .*sigma
	zzbins = minimum(zz) - zsig .* sigma : sigma : maximum(zz) + zsig .*sigma

	heigth  = 1000
	xxcontents = heigth * ones(Base.size(xx))

	coors = ndim == 2 ? (xx, yy) : (xx, yy, zz)
	edges = ndim == 2 ? (xxbins, yybins) : (xxbins, yybins, zzbins)
	hh  = SB.fit(SB.Histogram, coors, SB.weights(xxcontents), edges)

	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]

	factor = 10.0
	sigma_kernel = (sigma * factor) .* ones(ndim)
	weights_ = IF.imfilter(hh.weights, IF.Kernel.gaussian(sigma_kernel))

	cells    = findall(x -> x .> threshold, weights_)
	coors    = [[centers[i][cell[i]] for cell in cells] for i in 1:1:ndim]
	contents = [hh.weights[cell] for cell in cells]

	contents = [weights_[cell] for cell in cells]

	return (coors = coors, cells = cells, contents = contents, edges = edges)

	#img0 = jc.histo_to_img(h2)
	#sigma_kernel = 2 .* ones(ndin)
	#contents = IF.imfilter(hh.weights, IF.Kernel.gaussian(sigma_kernel))

	#h0 = SB.zero(hh)

	return hh
	#img = jc.Img(img0.coors, contents, img0.edges)
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
md"""

Notthing

"""

# ╔═╡ 6c8bf138-8fec-4c69-b4dd-4284faddeed0
begin
nndim = 3
img = line(ndim = nndim, threshold = 6.)
#sc2 = scatter3d(idf.x, idf.y, idf.z, marker_z = idf.energy,
#	markersize = vscale(idf.energy),
#	markertype = "circle", label = false, alpha = 0.4, c = :inferno,
#	xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")
#scatter3d!(imc.x, imc.y, imc.z, color = "white", alpha = 0.4, markersize = 0.5)
#plot(sc2, size = (700, 600))

#plot(img)
#typeof(img)
end

# ╔═╡ e2a822e9-10c2-481e-a770-f958f376a674
length(img.coors[1] )

# ╔═╡ f5674874-b238-4d0b-8df1-4e55e53fa446
histogram(img.contents, nbins = 100, title = "energy")

# ╔═╡ 009566ce-b9b5-4895-8565-6fb361845484
begin
theme(:dark)
scatter(img.coors..., zcolor = img.contents, alpha = 0.4, title = "energy")
end

# ╔═╡ fdcccc29-bf00-4f5f-aa85-ed86be2d4400
begin
if nndim == 2
	histogram2d(img.coors..., weights = img.contents, nbins = img.edges, title = "energy")
end
end

# ╔═╡ 5a1832c1-33ff-45dc-8f47-212179dbe862
md"""

Clouds 3D code

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
    #xborders, xneigh = _neighbour_node((xs, ys), xnodes, edges, steps, m)

    # output
    return (edges = edges,  coors = coors, contents = contents[cells],
			cells = cells,
            grad = grad[cells], igrad = igrad[cells],
            lap = lap[cells], curves = curves,
            curmax = curmax[cells], icurmax = icurmax[cells],
            curmin = curmin[cells], icurmin = icurmin[cells],
            nodes  = xnodes)
            #nborders = xborders, neighbour_node = xneigh)

end
end # begin

# ╔═╡ f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
begin
steps_ = [edge[2]-edge[1] for edge in img.edges]
xcl   = clouds(img.coors, img.contents, steps_)
end;

# ╔═╡ f00b443d-c5db-49cd-8df7-160df7b45f61
xcl.nodes

# ╔═╡ 13ac9fdf-46d0-4940-80e3-8619f0609108
md"""

## Plots
"""

# ╔═╡ 6fdb8681-c4a1-4e00-910c-8680fc65ff63
histogram(xcl.grad, nbins = 100, title = "gradient")

# ╔═╡ 702a11da-09f4-4159-ba64-cca34a81caa4
histogram(xcl.igrad, nbins = 50, title = "i-gradient")

# ╔═╡ ea72b2b1-3849-441e-832d-ae27bff5905a
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.grad, alpha = 0.1, title = "gradient")
end

# ╔═╡ 9311a861-4c21-4772-a2ea-3e9a7e330981
begin
if nndim == 2
	histogram2d(xcl.coors..., weights = xcl.grad, nbins = xcl.edges, title = "gradient")
end
end

# ╔═╡ 536dee15-9660-4887-b472-3f64da2a381b
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.igrad, alpha = 0.1, title = "gradient")
end

# ╔═╡ 327e8996-2bc0-4b85-9b57-e71ad8bc90d0
begin
if nndim == 2
	histogram2d(xcl.coors..., weights = xcl.igrad, nbins = xcl.edges, title = "gradient")
end
end

# ╔═╡ a66a215c-9598-47e1-876e-d8b670c2e2e4
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.lap, alpha = 0.1, title = "laplacian")
end

# ╔═╡ db1c0d58-1f33-4aa3-9be2-f4b0643e62a2
begin
if nndim == 2
	histogram2d(xcl.coors..., weights = xcl.lap, nbins = xcl.edges, title = "laplacian")
end
end

# ╔═╡ 60b83c6d-7ff5-4ceb-8a71-b51be8f2c157
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.curmin, alpha = 0.1, title = "cur min")
end

# ╔═╡ b14bfefe-251d-42ab-ae7e-67f8ca9ad141
begin
if nndim == 2
	histogram2d(xcl.coors..., weights = xcl.curmin, nbins = xcl.edges, title = "laplacian")
end
end

# ╔═╡ 68d7c4e9-3b97-45e8-aefc-a038a97616ce
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.curmax, alpha = 0.1, title = "cur min")
end

# ╔═╡ 3811c9f9-1369-4393-b28d-e6a5c2776e24
begin
if nndim == 2
	histogram2d(xcl.coors..., weights = xcl.curmax, nbins = xcl.edges, title = "laplacian")
end
end

# ╔═╡ be2f8a15-6c8f-4e9e-8983-1c04421ab0c7
begin
theme(:dark)
scatter(xcl.coors..., zcolor = xcl.nodes, alpha = 0.1, title = "nodes")
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
# ╠═e2a822e9-10c2-481e-a770-f958f376a674
# ╠═f5674874-b238-4d0b-8df1-4e55e53fa446
# ╠═009566ce-b9b5-4895-8565-6fb361845484
# ╠═fdcccc29-bf00-4f5f-aa85-ed86be2d4400
# ╠═5a1832c1-33ff-45dc-8f47-212179dbe862
# ╠═9c413a7e-3cbc-4676-b846-37db33f0c2f0
# ╠═1367811b-355c-446d-9d58-d07a71dfdb23
# ╠═23c6a3f5-443f-4199-87c8-68b93e495107
# ╠═f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
# ╠═f00b443d-c5db-49cd-8df7-160df7b45f61
# ╠═13ac9fdf-46d0-4940-80e3-8619f0609108
# ╠═6fdb8681-c4a1-4e00-910c-8680fc65ff63
# ╠═702a11da-09f4-4159-ba64-cca34a81caa4
# ╠═ea72b2b1-3849-441e-832d-ae27bff5905a
# ╠═9311a861-4c21-4772-a2ea-3e9a7e330981
# ╠═536dee15-9660-4887-b472-3f64da2a381b
# ╠═327e8996-2bc0-4b85-9b57-e71ad8bc90d0
# ╠═a66a215c-9598-47e1-876e-d8b670c2e2e4
# ╠═db1c0d58-1f33-4aa3-9be2-f4b0643e62a2
# ╠═60b83c6d-7ff5-4ceb-8a71-b51be8f2c157
# ╠═b14bfefe-251d-42ab-ae7e-67f8ca9ad141
# ╠═68d7c4e9-3b97-45e8-aefc-a038a97616ce
# ╠═3811c9f9-1369-4393-b28d-e6a5c2776e24
# ╠═be2f8a15-6c8f-4e9e-8983-1c04421ab0c7
# ╠═068e9533-1e4a-40be-83db-617a17935b0c
# ╠═654fb60f-2349-4b22-934e-dfb43080f5ec
