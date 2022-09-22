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

# ╔═╡ a57cdb41-c388-4976-bec8-ec0650fb139c
import jclouds as jc

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

# ╔═╡ f5fb4efe-3d08-4d75-8855-670d0fe92f25
begin
ximg = jc.cc_braid()
end

# ╔═╡ 6faf5d14-f3f2-4549-bc1b-7d46771eb1da
begin
xm, ym, zz   = vec(ximg.coors[1]), vec(ximg.coors[2]), vec(ximg.contents)
xbins, ybins = ximg.edges[1], ximg.edges[2]	
h1 = fit(Histogram, (xm, ym), weights(zz), (xbins, ybins))
plot(h1)
end

# ╔═╡ 8f49762d-9e98-4fed-b55c-17d9425d0fff
histogram2d(vec(xm), vec(ym), weights = vec(zz), nbins = (xbins, ybins))

# ╔═╡ 62867032-c903-4b5b-b6db-51c02ddc6e6d
begin
ximg2 = jc.histo_to_img(h1)
x2, y2, z2 = vec(ximg2.coors[1]), vec(ximg2.coors[2]), vec(ximg2.contents) 
histogram2d(x2, y2, weights = z2, bins = ximg2.edges)
#histogram2d(ximg2.coors..., weights = ximg2.contents, bins = ximg2.edges)
end

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

# ╔═╡ 441a1b38-7ff8-456c-8511-98f57828b26b
begin
	steps_   = [1., 1.]
	coors_   = (vec(ximg.coors[1]), vec(ximg.coors[2]))
	weights_ = vec(ximg.contents)
	mm       = jc.cc_moves2d()
	cl       = jc.cc_clouds(coors_, weights_, steps_)
end;

# ╔═╡ 7e1adeb7-2714-4205-b895-62cbe0477008
keys(cl)

# ╔═╡ 7856646f-20f1-47b1-986b-0702e2f53305
md"""
### Control Histograms
"""

# ╔═╡ bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
begin
xs = vec(cl.coors[1])
ys = vec(cl.coors[2])
end;

# ╔═╡ b774e74d-0843-4f5c-98b9-9103bd0a5657
histogram2d(xs, ys, weights = vec(cl.contents), 
		bins = cl.edges, alpha = 0.8, title = "contents")
	

# ╔═╡ a684c917-3750-4d80-a86e-7cc3fd7b9d02
begin
	histogram2d(xs, ys, weights = vec(cl.grad), 
		bins = cl.edges, alpha = 0.8, title = "gradient")
	quiver!(xs, ys, quiver = jc.cc_quiver(cl.igrad, mm, steps_))
end

# ╔═╡ d8d579a5-aa12-4128-b2ff-252d165cd2a6
begin
	histogram2d(xs, ys, weights = vec(cl.igrad), 
		bins = cl.edges, alpha = 0.8, title = "gradient direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.igrad, mm, steps_))
end

# ╔═╡ 1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
histogram2d(xs, ys, weights = vec(cl.lap), 
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
chis2 = [histogram2d(xs, ys, weights = vec(cl.curves[m]), 
		bins = cl.edges, title = m) for m in mm.moves]
plot(chis2..., layout = (3, 3))
end

# ╔═╡ c7dbc705-c905-4fc8-b1a5-616345d029b8
begin
	histogram2d(xs, ys, weights = vec(cl.curmax), 
		bins = cl.edges, alpha = 0.5, title = "curmax")
	quiver!(xs, ys, quiver = jc.cc_quiver(cl.icurmax, mm, steps_))
end

# ╔═╡ 62d11578-0d90-472d-8898-83e21c53d621
begin
	histogram2d(xs, ys, weights = vec(cl.icurmax), 
		bins = cl.edges, alpha = 0.5, title = "curmax direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmax, mm, steps_))
end

# ╔═╡ 1605c6b4-f674-4209-ae46-f7ac4813693d
begin
	histogram2d(xs, ys, weights = vec(cl.curmin), 
		bins = cl.edges, alpha = 0.5, title = "curmin direction")
	quiver!(xs, ys, quiver = quiver(cl.icurmin, mm, steps_))
end

# ╔═╡ f4155f75-af7f-4b79-b846-3bdc601d8767
begin
	histogram2d(xs, ys, weights = vec(cl.icurmin), 
		bins = cl.edges, alpha = 0.5, title = "curmin direction")
	#quiver!(vec(cl.x), vec(cl.y), quiver = quiver(cl.icurmin, mm, steps_))
end

# ╔═╡ 654bce47-b51d-4c9a-81b3-b295f4bb055a
begin
	histogram2d(xs, ys, weights = vec(cl.nodes), 
		bins = cl.edges, alpha = 1., title = "nodes")
end

# ╔═╡ 503b611f-d2f9-4b94-92c5-0b005505d5bf
begin
	histogram2d(xs, ys, weights = vec(cl.nborders), 
		bins = cl.edges, alpha = 1., title = "number of borders")
end

# ╔═╡ 27f27825-3f01-4c63-a119-4267ef69b11c
begin
nhis2 = [histogram2d(xs, ys, weights = vec(cl.neighbour_node[i]), 
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
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╠═3922eba2-f322-4b06-b9e0-83bc723d7930
# ╠═3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╠═ceaff28e-ea37-4a92-a812-210fd2b91fad
# ╠═94bb1cff-3f57-4af6-b177-ff732cc78429
# ╠═7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╠═40a8cad8-0a7b-41c5-9807-a5143c30a76b
# ╠═f5fb4efe-3d08-4d75-8855-670d0fe92f25
# ╠═6faf5d14-f3f2-4549-bc1b-7d46771eb1da
# ╠═8f49762d-9e98-4fed-b55c-17d9425d0fff
# ╠═62867032-c903-4b5b-b6db-51c02ddc6e6d
# ╠═bc0a9ad4-62f4-4ef9-909a-4a2f68584b6f
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
