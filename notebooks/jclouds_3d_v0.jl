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
#using LinearAlgebra
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
img = line(ndim = 3, threshold = 6.)
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
histogram(img.contents, nbins = 100)

# ╔═╡ 009566ce-b9b5-4895-8565-6fb361845484
begin
theme(:dark)
scatter(img.coors..., zcolor = img.contents, alpha = 0.4)
end

# ╔═╡ c76ba535-01c2-4403-a605-6164a5b6aa89
histogram2d(img.coors..., weights = img.contents, nbins = img.edges)

# ╔═╡ 7a0e6e69-9f9b-427c-8077-172fe9c8af7a
histogram(img.contents, nbins = 100)

# ╔═╡ 58a19a9e-1955-4a0e-8de8-a3e5d00d2340
begin
ndim = length(img.coors)
print(ndim)
#if ndim == 2
#theme(:dark)
#sc = scatter(img.coors..., zcolor = img.contents)
#plot(sc)
#elseif ndim == 3
#	sc = scatter3d(img.coors..., marker_z = img.contents,
#		markersize = vscale(img.contents),
#		markertype = "circle", label = false, alpha = 0.4, c = :inferno,
#		xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")
#end
end #begin

# ╔═╡ 9b96e42d-3c15-4863-b810-de19f2b89ce3
#sc2 = scatter3d(img.coors..., marker_z = img.contents,
#	markersize = vscale(img.contents),
#	markertype = "circle", label = false, alpha = 0.4, c = :inferno,
#	xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")

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

# ╔═╡ 602c48e9-9ba6-4cb2-b9f9-f9a5ed340622
img.edges

# ╔═╡ 17d4089d-ed31-4a1a-a49a-669b50b0f7b1
mesh(img.edges)

# ╔═╡ 654fb60f-2349-4b22-934e-dfb43080f5ec
begin
function _mesh3d(x, y, z)
	xv = getindex.(Iterators.product(x, y, z), 1)  # first.(Iterators.product(x, y, z), 1) is also ok
	yv = getindex.(Iterators.product(x, y, z), 2)
	zv = getindex.(Iterators.product(x, y, z), 3)
	return xv, yv, zv
end
end #begin

# ╔═╡ 458280ed-fc18-48e5-844e-f6e21dc1ad19
mask = img.weights .> 0

# ╔═╡ b26fca8e-6593-496d-9ff2-3ecbeea3f69f
indices = CartesianIndices(img.weights .> 0)

# ╔═╡ 03aa6f7f-c30f-4ac2-a51f-f559cde10247
begin
cells = findall(x -> x.>0, img.weights)
_xs   = [img.edges[1][cell[1]] for cell in cells]
end

# ╔═╡ bea99be8-5f32-4e08-bdc9-ef1e8b0746eb
cell0[1]

# ╔═╡ 8d7d18f8-5c8a-4762-bb07-90d6d84cf07f
plot(img)

# ╔═╡ bd311cb5-6e59-4958-9755-7468f1faf0fd
md"""

## Dev
"""

# ╔═╡ 1d8af6cd-c48b-4cbe-a040-70eb4dfbc0b4
begin
histogram2d(vec(img.coors[1]), vec(img.coors[2]), weights = vec(img.contents), nbins = img.edges)
end

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

# ╔═╡ bc0a9ad4-62f4-4ef9-909a-4a2f68584b6f
md"""

## Clouds

Cloud creates an internal image using as input the coordinates (x, y) the weights and the steps (xstep, ystep)

It computes the discrete local gradient, it direction, the laplacian and the maximun and minimum curvature and its direction. Only weights > 0 are consider.

It also groups the cells into nodes, and connect the nodes if both have connecting cells.

It returns a NamedTuple with all the previous information.

"""

# ╔═╡ 1cce8566-ad73-4b6c-abb9-47d95c17e412
begin
typeof(img)
end

# ╔═╡ 441a1b38-7ff8-456c-8511-98f57828b26b
begin
	steps    = keyimg == :line ? [0.5, 0.5] : [1.0, 1.0]
	coors    = (vec(img.coors[1]), vec(img.coors[2]))
	contents = vec(img.contents)
	mm       = jc.moves2d()
	cl       = jc.clouds(coors, contents, steps)
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
	quiver!(xs, ys, quiver = jc.quiver(cl.igrad, mm, steps))
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
	quiver!(xs, ys, quiver = jc.quiver(cl.icurmax, mm, steps_))
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
	#quiver!(xs, ys, quiver = jc.quiver(cl.icurmin, mm, steps))
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

function cloud_graph(nodes, neighs)
	dus = jc._edges(nodes, neighs)

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
# ╠═c76ba535-01c2-4403-a605-6164a5b6aa89
# ╠═7a0e6e69-9f9b-427c-8077-172fe9c8af7a
# ╠═58a19a9e-1955-4a0e-8de8-a3e5d00d2340
# ╠═9b96e42d-3c15-4863-b810-de19f2b89ce3
# ╠═068e9533-1e4a-40be-83db-617a17935b0c
# ╠═602c48e9-9ba6-4cb2-b9f9-f9a5ed340622
# ╠═17d4089d-ed31-4a1a-a49a-669b50b0f7b1
# ╠═654fb60f-2349-4b22-934e-dfb43080f5ec
# ╠═458280ed-fc18-48e5-844e-f6e21dc1ad19
# ╠═b26fca8e-6593-496d-9ff2-3ecbeea3f69f
# ╠═03aa6f7f-c30f-4ac2-a51f-f559cde10247
# ╠═bea99be8-5f32-4e08-bdc9-ef1e8b0746eb
# ╠═8d7d18f8-5c8a-4762-bb07-90d6d84cf07f
# ╠═bd311cb5-6e59-4958-9755-7468f1faf0fd
# ╠═1d8af6cd-c48b-4cbe-a040-70eb4dfbc0b4
# ╠═40a8cad8-0a7b-41c5-9807-a5143c30a76b
# ╠═bc0a9ad4-62f4-4ef9-909a-4a2f68584b6f
# ╠═1cce8566-ad73-4b6c-abb9-47d95c17e412
# ╠═441a1b38-7ff8-456c-8511-98f57828b26b
# ╠═7e1adeb7-2714-4205-b895-62cbe0477008
# ╠═7856646f-20f1-47b1-986b-0702e2f53305
# ╠═bd1711b7-b6c3-4a97-a54e-a7c1e268ef92
# ╠═b774e74d-0843-4f5c-98b9-9103bd0a5657
# ╠═a684c917-3750-4d80-a86e-7cc3fd7b9d02
# ╠═d8d579a5-aa12-4128-b2ff-252d165cd2a6
# ╠═1cf114b7-3f8e-4f4c-bdcd-81e0b6ddee74
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
