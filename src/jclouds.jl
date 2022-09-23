module jclouds

using Base
using StatsBase
using LinearAlgebra

export Img, Moves2D
export braid, moves2d
export mesh, quiver
export clouds
export _edges

#-------
# Imgs
#-------

struct Img
    coors
    contents
    edges
end

"""
Create a braid Img (trenza) in 2D
return imgs: (coors, contents, edges)
"""
function braid(nn = 10)
	nn = 10
	zz = 0:nn-1
	zz = vcat(zz, reverse(zz))
	zz = ones(nn)' .* zz
	xx = 0.5 .+ 0:2*nn
	yy = 0.5 .+ 0:nn
	xm = ones(nn)' .* xx
	ym = yy'       .* ones(2*nn)
	xbins = 0:1:2*nn
	ybins = 0:1:nn
    return Img((xm, ym), zz, (xbins, ybins))
	#return (coors = (xm, ym), contents = zz, edges = (xbins, ybins))
end

"""
convert a 2D histogram into img
"""
function histo_to_img(h::Histogram)
	nx, ny  = Base.size(h.weights)
	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in h.edges]
	xm = ones(ny)'   .* centers[1]
	ym = centers[2]' .* ones(nx)
	zz = deepcopy(h.weights)
	return Img((xm, ym), zz, h.edges)
end

#-----------
# Moves
#-----------


struct Moves2D
	moves
	imove0
	dmoves
	omoves
end

function moves2d()
	moves  = [[i, j] for i in -1:1:1 for j in -1:1:1]
	imove0 = [i for (i, move) in enumerate(moves) if move == [0, 0]][1]
	dmoves = Dict(1:9 .=> moves)
	omoves = Dict()
	imoves = deepcopy(moves)
	imoves = filter(move -> move != imoves[imove0], moves)
	for i in 1:9
		omoves[moves[i]] = [imove for imove in imoves
			if LinearAlgebra.dot(imove, moves[i]) == 0]
	end
	str_moves = Moves2D(moves, imove0, dmoves, omoves)
end


#-------
#  Clouds
#-------

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

#--- Nodes

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

#--- Edged

function _edges(nodes, neighs)
	iinodes = sort(unique(vec(nodes[nodes .>0])))
	dus     = Dict()
	for iii in iinodes
		imask = nodes .== iii
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
	end
	return dus
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

function clouds(coors, energy, steps)

    ndim  = length(coors)
    nsize = length(coors[1])

    # assert dimensions
    for i in 2:ndim
        @assert(length(coors[2]) == nsize)
    end
    @assert(length(energy) == nsize)
    @assert(length(steps)  == ndim)

    # define the extended edges
    edges = [minimum(x)-1.5*step:step:maximum(x)+1.5*step for (x, step) in zip(coors, steps)]

    # alias
    xx, yy         = coors
    xedges, yedges = edges

    # main histogram
    histo    = fit(Histogram, (xx, yy), weights(energy), (xedges, yedges))
    contents = deepcopy(histo.weights)

    # deltas
    m = moves2d()
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

    # output
    return (edges = (xedges, yedges),  coors = (xs, ys), contents = contents,
            grad = grad, igrad = igrad,
            lap = lap, curves = curves,
            curmax = curmax, icurmax = icurmax,
            curmin = curmin, icurmin = icurmin,
            nodes  = xnodes,
            nborders = xborders, neighbour_node = xneigh)

end


end # end of module
