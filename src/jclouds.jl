module jclouds

using Base
import StatsBase as SB
import LinearAlgebra as LA


#export mesh, quiver
export line, clouds


#-----------
# Moves
#-----------

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



#-----------------------------
#  Clouds internal functions
#-----------------------------

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
		if (i != m.imove0)
        	dd = curves[move]
        	mask1 = dd .> curmax
        	curmax[mask1] .= dd[mask1]
        	icurmax[mask1] .= i
        	mask2 = dd .< curmin
        	curmin[mask2] .= dd[mask2]
        	icurmin[mask2] .= i
		end
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


#--- Public
#
# function mesh(edges)
#     nx, ny  = length(edges[1])-1, length(edges[2])-1
#     centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]
#     xm = ones(ny)'   .* centers[1]
#     ym = centers[2]' .* ones(nx)
#     return xm, ym
# end
#

# function _mesh3d(x, y, z)
# 	xv = getindex.(Iterators.product(x, y, z), 1)  # first.(Iterators.product(x, y, z), 1) is also ok
# 	yv = getindex.(Iterators.product(x, y, z), 2)
# 	zv = getindex.(Iterators.product(x, y, z), 3)
# 	return xv, yv, zv
# end

#
# function quiver(idir, m, steps)
#     xstep, ystep = steps
#     uus = [m.dmoves[k] for k in vec(idir)]
#     us = [xi * 0.8 * xstep for (xi, yi) in uus]
#     vs = [yi * 0.8 * ystep for (xi, yi) in uus]
#     return us, vs
# end
#
# """
# linear scale a vector of values between a minimum and a maximum
#
# Parameters:
# 	var : Vector{Real}
# 	emin: Real
# 	emax: Real
#
# Return:
# 	vvar: Vector{Real}
#
# """
# function vscale(var, emin = 1, emax = 4)
# 	vvar = (var .- minimum(var)) ./(maximum(var) - minimum(var)) .*(emax-emin) .+ emin
# 	return vvar
# end

end # end of module
