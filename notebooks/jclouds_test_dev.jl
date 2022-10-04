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


Test clouds 3D in Julia


J.A. Hernado,

Santiago, October 2022

---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly();

# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""

## Generate Matrix

"""

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin
bndim = @bind ndim Select([2, 3])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select dimensions of the line $(bndim)

"""
end

# ╔═╡ 34820285-e673-4f45-9593-fc5cb409d3d1
begin

function get_box(vals = [1, 2, 4])

	@assert (length(vals) >= 3) & (length(vals) <= 4)
	ndim = length(vals) - 1
	
	a, b, c = vals[1], vals[2], vals[3]
	@assert (a <= b) & (b <= c)
	d = ndim == 3 ? vals[3] : 0
	
	cs = [a b a; b c b; a b a]
	
	bs = zeros(3, 3, 3)
	for k in 1:3
		for i in 1:3
			for j in 1:3
			offset = k == 2 ? d : 0
			bs[i, j, k] = cs[i, j] + offset
			end
		end
	end

	mat = ndim == 2 ? cs : bs

	fx = b-a
	fy = c-b
	fxy = (c-a)/sqrt(2.)

	c1, c2, c3 = max(fx, fxy), fy, 0.
	grad = [c1 c2 c1; c2 c3 c2; c1 c2 c1]

	lx = 2*fx + fxy 
	ly = -2*fx + fy  
	lz = -4*fxy - 4*fy
	lap = [lx ly lx; ly lz ly; lx ly lx]

	c3 = 2*min(fxy, fy)
	maxc = [c1 c2 c1; c2 c3 c2; c1 c2 c1]

	c1, c2, c3 = min(fx, fxy, 0.0), -2*fx, min(-2*fy, -2*fxy)
	minc = [c1  c2 c1; c2 c3 c2; c1 c2 c1]
	
	return mat, grad, lap, maxc, minc
	
end

end

# ╔═╡ 43006244-5f3c-4793-810c-a46189028a6f
mat, grad, lap, curmax, curmin = get_box([1, 2, 10]) 

# ╔═╡ b8d9d2b9-543c-44f9-bb69-b9c3c8687145
min(2, 3)

# ╔═╡ aec7f4d6-ec2d-4646-bdde-153d85005d45
md"""

## Test code

"""

# ╔═╡ 7e32f198-6473-4399-9dc4-af54ed451c33
begin
function _test_clouds(bs, threshold = 0.0)

	# set the input for clouds
	ndim   = length(size(bs)) 
	indices = CartesianIndices(bs)
	coors = [vec([index[i] for index in indices]) for i in 1:ndim]
	steps = ones(ndim)

	# call clouds
	tcl = jc.clouds(coors, vec(bs), steps)

	# prepare the local matrix to compute gradients, laplacian and curvatures
	nsize  = size(bs) .+ 2	
	aa    = zeros(Float64, nsize...)
	for index in indices
		kindex = CartesianIndex(Tuple(index) .+ 1)
		aa[kindex] = bs[index]
	end

	# the mask
	mask = aa  .>  threshold
	nmask = aa .<= threshold
	# set the moves
	mm = jc.moves(ndim)

	# compute the deltas
	function _delta(move)
		delta = zeros(Float64, nsize...)
		mod   = sqrt(sum(move .* move))
		mod   = mod == 0.0 ? 1 : mod
		for index in indices
			uindex = CartesianIndex(Tuple(index)  .+ 1)
			kindex = CartesianIndex(Tuple(Tuple(uindex) .+ move))
			dd = aa[kindex] > 0.0 ? (aa[kindex] - aa[uindex])/mod : 0.0
			delta[uindex] = dd
		end
		return delta
	end

	# compute deltas
	deltas = [_delta(move) for move in mm.moves]

	# compute gradient
	ugrad = zeros(Float64, nsize...)
	igrad = mm.i0 .* ones(Int, nsize...)
	for (i, delta) in enumerate(deltas)
		imask = delta .> ugrad 
		ugrad[imask] = delta[imask]
		igrad[imask] .= i
	end

	# test the gradient
	@assert sum(ugrad[mask] .== tcl.grad)  == sum(mask)
	@assert sum(igrad[mask] .== tcl.igrad) == sum(mask)
	print("Ok grad! \n")

	# compute and test laplacian
	lap = reduce(.+, deltas)
	@assert sum(lap[mask] .== tcl.lap) == sum(mask)
	print("Ok laplacian! \n")

	# compute and test curves and max, min curve
	curves = [reduce(.+, [deltas[i] for i in mm.iortho[s]]) for s in mm.isym]
	curmax = -1e6*ones(nsize...)
	curmin = +1e6*ones(nsize...)
	icurmax = mm.i0 .* ones(Int, nsize...)
	icurmin = mm.i0 .* ones(Int, nsize...)
	for (i, curve) in enumerate(curves)
		imove = mm.isym[i][1]
		imask = curve .> curmax 
		curmax[imask] .= curve[imask]
		icurmax[imask] .= imove
		imask = curve .< curmin
		curmin[imask] .= curve[imask]
		icurmin[imask] .= imove
	end
	icurmin[nmask] .= mm.i0
	icurmax[nmask] .= mm.i0

	@assert sum(curmin[mask] .== tcl.curmin) == sum(mask)
	@assert sum(curmax[mask] .== tcl.curmax) == sum(mask)
	@assert sum(icurmin[mask] .== tcl.icurmin) == sum(mask)
	@assert sum(icurmax[mask] .== tcl.icurmax) == sum(mask)
	print("OK curvmax/min icurvmax/min \n")
	

	#return icurmin[mask], tcl.icurmin

	return tcl
	#return curmax, curmin, icurmax, icurmin
	
end
	
end

# ╔═╡ b0107bba-91e3-4c66-a1a9-b5bed75b308c
tcl = _test_clouds(mat)

# ╔═╡ 2f98878b-75f5-4a94-9f1f-aa35b8843f25
begin
@assert sum(tcl.grad   .== vec(grad))   == 9
@assert sum(tcl.lap    .== vec(lap))    == 9
#@assert sum((tcl.curmax, 3, 3) .== curmax) == 9
#@assert sum((tcl.curmin, 3, 3) .== curmin) == 9
end

# ╔═╡ f0d31ed4-9a64-4fd1-8781-dd68879fc8eb
reshape(tcl.grad, 3, 3), grad

# ╔═╡ ca59c78f-f04f-438b-9b95-c8826dd1890f
reshape(tcl.lap, 3, 3), lap

# ╔═╡ 6c763593-a692-442e-936d-5a8bba3ef185
reshape(tcl.curmax, 3, 3), curmax

# ╔═╡ 4118b1cb-53b8-492d-a09e-ea5275182eb5
reshape(tcl.curmin, 3, 3), curmin

# ╔═╡ c6a33da6-446f-49bb-acd4-6b024538429e
function xx_test_clouds(cs)

	# prepare clouds inputs
	nx, ny = size(cs)
	is = vcat([[i for j in 1:ny] for i in 1:nx]...)
	js = repeat([j for j in 1:ny], nx)
	steps = [1.0, 1.0]

	# create clouds
	tcl = jc.clouds((is, js), vec(cs), steps)

	# compute extended array 
	nx, ny = size(cs)
	aa = zeros(Float64, nx + 2, ny + 2)
	for i in 1:nx
		for j in 1:ny
			aa[i+1, j+1] = cs[i, j] 
		end
	end
	asize = size(aa)

	#moves 
	moves = [[i, j] for i in -1:1:1 for j in -1:1:1]

	# compute deltas
	function _delta(aa, move)
		ai = zeros(Float64, asize...)
		for i in 2:4
			for j in 2:4
				di, dj = move[1], move[2]
				dmove = sqrt(di*di + dj*dj) 
				dmove = dmove == 0. ? 1. : dmove
				dd = aa[i+di, j + dj] >= 0.0 ? (aa[i + di, j + dj] - aa[i, j])/dmove : 0
				ai[i, j] = dd
			end
		end
		return ai
	end
	deltas = [_delta(aa, move) for move in moves]

	# compute gradient
	ugrad = zeros(asize...)
	for delta in deltas
		mask = delta .> ugrad
		ugrad[mask] = delta[mask]
	end
	ugrad

	# test gradient
	mask = aa .> 0.0
	@assert sum(ugrad[mask] .== tcl.grad) == sum(mask)

	# compute and test laplacian
	lap = reduce(.+, deltas)
	@assert sum(lap[mask] .== tcl.lap) == sum(mask)

	# compute and test curves and max, min curve
	curves = [deltas[i] .+ deltas[j] for (i, mi) in enumerate(moves) for (j, mj) in enumerate(moves) if (sum(mi .+ mj) == 0.0) & (sum(mi .* mj) != 0.0) & ( j > i)]
	curmax = -1e6*ones(asize...)
	curmin = +1e6*ones(asize...)
	for curve in curves
		imask = curve .> curmax
		curmax[imask] .= curve[imask]
		imask = curve .< curmin
		curmin[imask] .= curve[imask]
	end
	@assert sum(curmin[mask] .== tcl.curmin) == sum(mask)
	@assert sum(curmax[mask] .== tcl.curmax) == sum(mask)

	return true
end

# ╔═╡ adf969fd-1555-4099-a197-50281e66597e
begin
function _test_clouds_new(bs, threshold = 0.0)

	# set the input for clouds
	ndim   = length(size(bs)) 
	indices = CartesianIndices(bs)
	coors = [vec([index[i] for index in indices]) for i in 1:ndim]
	steps = ones(ndim)

	# call clouds
	tcl = jc.clouds(coors, vec(bs), steps)
	

	# prepare the local matrix to compute gradients, laplacian and curvatures
	nsize  = size(bs) .+ 2	
	aa    = zeros(Float64, nsize...)
	for index in indices
		kindex = CartesianIndex(Tuple(index) .+ 1)
		aa[kindex] = bs[index]
	end

	# the mask
	mask = aa  .>  threshold
	nmask = aa .<= threshold
	aa[nmask]  .=  threshold
	
	# set the moves
	mm = _moves(ndim)

	# compute the deltas
	function _delta(move)
		delta = zeros(Float64, nsize...)
		mod   = sqrt(sum(move .* move))
		mod   = mod == 0.0 ? 1 : mod
		for index in indices
			uindex = CartesianIndex(Tuple(index)  .+ 1)
			kindex = CartesianIndex(Tuple(Tuple(uindex) .+ move))
			#dd = mask[kindex] ? (aa[kindex] - aa[uindex])/mod : 0.0
			dd = true ? (aa[kindex] - aa[uindex])/mod : 0.0
			delta[uindex] = dd
		end
		return delta
	end

	# compute deltas
	deltas = [_delta(move) for move in mm.moves]
	return deltas

	# compute gradient
	ugrad = zeros(Float64, nsize...)
	igrad = mm.i0 .* ones(Int, nsize...)
	for (i, delta) in enumerate(deltas)
		imask = delta .> ugrad 
		ugrad[imask] = delta[imask]
		igrad[imask] .= i
	end

	# test the gradient
	@assert sum(ugrad[mask] .== tcl.grad)  == sum(mask)
	@assert sum(igrad[mask] .== tcl.igrad) == sum(mask)
	print("Ok grad! \n")

	# compute and test laplacian
	lap = reduce(.+, deltas)
	@assert sum(lap[mask] .== tcl.lap) == sum(mask)
	print("Ok laplacian! \n")

		# compute and test curves and max, min curve
	curves = [deltas[i] .+ deltas[j] for (i, j) in mm.isym]
	curmax = -1e6*ones(nsize...)
	curmin = +1e6*ones(nsize...)
	icurmax = mm.i0 .* ones(Int, nsize...)
	icurmin = mm.i0 .* ones(Int, nsize...)
	for (i, curve) in enumerate(curves)
		imove = mm.isym[i][1]
		imask = curve .> curmax 
		curmax[imask] .= curve[imask]
		icurmax[imask] .= imove
		imask = curve .< curmin
		curmin[imask] .= curve[imask]
		icurmin[imask] .= imove
	end
	icurmin[nmask] .= mm.i0
	icurmax[nmask] .= mm.i0

	@assert sum(curmin[mask] .== tcl.curmin) == sum(mask)
	@assert sum(curmax[mask] .== tcl.curmax) == sum(mask)
	#@assert sum(icurmin[mask] .== tcl.icurmin) == sum(mask)
	#@assert sum(icurmax[mask] .== tcl.icurmax) == sum(mask)


	return icurmin[mask], tcl.icurmin
	
	#return curmax, curmin, icurmax, icurmin
	
end
	
end

# ╔═╡ Cell order:
# ╠═5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╠═a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─3922eba2-f322-4b06-b9e0-83bc723d7930
# ╟─3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╟─7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╟─e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╠═34820285-e673-4f45-9593-fc5cb409d3d1
# ╠═43006244-5f3c-4793-810c-a46189028a6f
# ╠═b8d9d2b9-543c-44f9-bb69-b9c3c8687145
# ╟─aec7f4d6-ec2d-4646-bdde-153d85005d45
# ╠═7e32f198-6473-4399-9dc4-af54ed451c33
# ╠═b0107bba-91e3-4c66-a1a9-b5bed75b308c
# ╠═2f98878b-75f5-4a94-9f1f-aa35b8843f25
# ╠═f0d31ed4-9a64-4fd1-8781-dd68879fc8eb
# ╠═ca59c78f-f04f-438b-9b95-c8826dd1890f
# ╠═6c763593-a692-442e-936d-5a8bba3ef185
# ╠═4118b1cb-53b8-492d-a09e-ea5275182eb5
# ╠═c6a33da6-446f-49bb-acd4-6b024538429e
# ╠═adf969fd-1555-4099-a197-50281e66597e
