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

## Generate Image

Produces a 3x3 box where to check the values of grad, lap, and maxcur and mincur

"""

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin
bndim = @bind nndim Select([2, 3])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select dimensions of the line $(bndim)

"""
end

# ╔═╡ 9e00d479-f411-4dfd-b51f-95dc9efeecd1
begin
cs = [1 2 1; 2 4 2; 1 2 1]
is = [ 1 1 1; 2 2 2; 3 3 3]
js = [1 2 3; 1 2 3; 1 2 3]
end

# ╔═╡ a6939579-e310-4ed9-a76d-85d09cba4ce1
begin
aa = zeros(5, 5)
for i in 2:4
	for j in 2:4
		aa[i, j] = cs[i-1, j-1]
	end
end
aa
end

# ╔═╡ df582445-6873-4829-9cce-82c1fb737b4b
length(size(aa))

# ╔═╡ c6a33da6-446f-49bb-acd4-6b024538429e
function _test_clouds(cs)

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
				ai[i, j] = (aa[i + di, j + dj] - aa[i, j])/dmove
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

# ╔═╡ 988498e0-6bc9-4c4d-9467-a941a848a157
begin
ok = _test_clouds(cs)

md"""

**Clouds has passed the 2d test? $(ok)**

"""
end

# ╔═╡ 938208b9-9a88-455a-a093-0e9ff4900fe9
begin
bs = zeros(3, 3, 3)
for k in 1:3
	for i in 1:3
		for j in 1:3
		offset = k == 2 ? 2 : 0
		bs[i, j, k] = 2.0 * cs[i, j] + offset
		end
	end
end
bs
end

# ╔═╡ 2c0570a4-c5d0-4a96-9fd9-9618cc1da7f3
begin
function _moves(ndim)
	moves = ndim == 2 ? [[i, j] for i in -1:1:1 for j in -1:1:1] : [[i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1]
	move0 = ndim == 2 ? [0, 0] : [0, 0, 0]
	kmove0 = [i for (i, move) in enumerate(moves) if (move == move0)][1]
	smoves = [(i, j) for (i, movei) in enumerate(moves) for (j, movej) in enumerate(moves) if (movei == -1 .* movej ) &  (i > j)]
	omoves = Dict()
	for ii in smoves
		movei = moves[ii[1]]
		omoves[ii] = [j for (j, movej) in enumerate(moves) if ((sum(movej .* movei) == 0.0) & (sum(movej .* movej) != 0.0))]
	end
	return (moves = moves, i0 = kmove0, isym = smoves, iorto = omoves)
end
end

# ╔═╡ fcb0a412-ccd6-4406-98f5-cd60ae9f7845
mm = _moves(3)

# ╔═╡ 2f331294-fa21-431a-928c-0b4d302726ec
mm.moves[mm.i0]

# ╔═╡ c70469df-5f59-4ea2-94e0-a247966f0f13
mm.isym

# ╔═╡ 792af853-8220-47e8-8a3f-865decbbeafe
begin
for ii in keys(mm.iorto)
	print([mm.moves[ii[j]] for j in 1:2], " = , ", [mm.moves[kk] for kk in mm.iorto[ii]], '\n')
end
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

# ╔═╡ b8476bc1-8ca1-4539-bdbb-71ddf3fa378c
cs, bs

# ╔═╡ 4efe53aa-83b8-46fa-8418-20023d3df966
_test_clouds_new(cs, 0.)

# ╔═╡ fc36f3ff-0f90-428b-887c-4761706a3dca
m2 = _moves(2)

# ╔═╡ b4971b5e-a614-4498-a177-ea218ae10933
mm.moves[1]

# ╔═╡ 683618f0-6df4-414c-a1d4-6bc27911cc34
bs

# ╔═╡ eb1c838c-42f7-47d9-98d5-9028f834d7c1
ones(3)

# ╔═╡ 3638a83e-f986-4752-96e4-a85725387e56
[i for (i, move) in enumerate(xmoves) if (move == [0, 0, 0])][1]

# ╔═╡ 7349abcb-1462-43e2-98f5-e4e84d114748
xmoves[24]

# ╔═╡ 17bdc350-b372-40c7-a4b1-6bf5720e03db
6.0/sqrt(2)

# ╔═╡ fc10da1f-8e07-486b-adc5-db1dc4304448
xmoves[15]

# ╔═╡ d84f1dcc-ba31-4a03-ba60-dcaea8b24938
kmoves

# ╔═╡ f7035d4a-01aa-4e58-9761-e3e2e346ccd0
bs

# ╔═╡ 6b14de6e-b40e-4017-94dd-2d34cfa66000
xmoves[5]

# ╔═╡ d8626b01-5007-45f2-b47e-b955fc8858ee
CartesianIndices(bs) .+ xmoves[1]

# ╔═╡ abda946a-f1b3-4306-9f30-d710fe5e80bd
sum(xmoves[1] .* xmoves[2])

# ╔═╡ 696e4101-c218-4988-b2f5-0d617f386283
xmoves

# ╔═╡ 04b6be28-ed7a-49c6-82c0-1666645e3f98
xmoves[1] .* xmoves[2]

# ╔═╡ 3f03a5ef-ecf4-48dc-bd56-71d365f5e553
xaa

# ╔═╡ bbc9921b-8429-40f8-bc02-93b294f3f3d3
bs

# ╔═╡ a1300bcc-ff12-4614-91cf-efbc4791f71a
bs

# ╔═╡ 3862afd3-fb51-482e-bbe6-ce8d20ec895b
begin
i1  = xindices[1] 
ii1 = CartesianIndex(Tuple(i1) .+ 1)
i1, ii1, bs[i1], bs[ii1]
end

# ╔═╡ ef73a5ad-e7dc-4bfc-8c06-5680b3f28c85
bs[ii]

# ╔═╡ b22dec95-bc58-44fe-b29e-af335c54e710
xndim

# ╔═╡ 6f5acf59-2a84-4b59-8010-fbd07b769429
function _test_clouds_3d(cs)

	# prepare clouds inputs
	ndim   = length(size(cs))
	nsize  = size(cs) .+ 2
	aa    = zeros(Float64, nsize...)
	indices = CartesianIndices(aa)
	coors = [[index[i] for index in indices] for i in 1:ndim]
	
	ndim   = length(cs)
	asize  = size(cs) .+ 2
	aa     = zeros(Float64, nsize...)
	
	coors  = zeros()
	is = vcat([[i for j in 1:ny] for i in 1:nx]...)
	js = repeat([j for j in 1:ny], nx)
	steps = [1.0, 1.0]

bcoors = [[index[i] for index in indices] for i in 1:3]
	
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
				ai[i, j] = (aa[i + di, j + dj] - aa[i, j])/dmove
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

# ╔═╡ a0c6aa94-ab32-4d92-a4f9-e17bc7682297
[i*j for i in 1:3 for j in 1:3 if (i == 2) & (j ==2)]

# ╔═╡ f0fe6965-fdbe-4f94-aac3-a41f3cdd63f0
moves = [[i, j] for i in -1:1:1 for j in -1:1:1]

# ╔═╡ 76f36865-51d6-4eb5-a242-677aa57427e0
begin
deltas = Dict()
function _delta(aa, move)
	ai = zeros(Float64, 5, 5)
	for i in 2:4
		for j in 2:4
			di, dj = move[1], move[2]
			dmove = sqrt(di*di + dj*dj) 
			dmove = dmove == 0. ? 1. : dmove
			ai[i, j] = (aa[i + di, j + dj] - aa[i, j])/dmove
		end
	end
	return ai
end
deltas = [_delta(aa, move) for move in moves]
end

# ╔═╡ c66c5247-09de-426b-a9d3-ffb2f69c5edf
_delta(xmoves[5])

# ╔═╡ d0f3e7c4-d72f-4822-b683-60ba5767278b
begin
curves = [ deltas[i] .+ deltas[j] for (i, mi) in enumerate(moves) for (j, mj) in enumerate(moves) if sum(mi .+ mj) == 0.0]
maxcurve = -1e6*ones(size(deltas[1])...)
for curve in curves
	mask = curve .> maxcurve
	maxcurve[mask] .= curve[mask]
end
maxcurve
end # begin

# ╔═╡ 62e67f7b-bc7a-4b2c-9800-7f89e7ef8d81
begin
ugrad = zeros(5, 5)
for delta in deltas
	mask = delta .> ugrad
	ugrad[mask] = delta[mask]
end
ugrad
end

# ╔═╡ b01800e5-18b8-4e3a-b1c8-34ada12d9854
reshape(ugrad[aa .> 0], 3, 3)

# ╔═╡ 36d0c79d-94b1-4454-b605-50791e51d184
ulap = reduce(.+, deltas)

# ╔═╡ 13bc2859-55b0-4670-9e42-d42e65d989f0


# ╔═╡ 3526ed61-2dc9-49ca-83eb-4feafd29e0ae
deltas

# ╔═╡ a796bc34-b937-4f71-a12b-baa5cbd266b9
begin
ai = zeros(Float64, 5, 5)
ai[2, 2] = 1
ai
end

# ╔═╡ 719e52ae-68fe-43c9-877e-1d0e9043e4b7
Vector{Array{Float64, (2, 2)}}

# ╔═╡ 3bbc7b6f-bfa2-4479-bfbc-93dbd695bdbb


# ╔═╡ 34923e48-1392-4c29-8c77-f1c58c89a3e2
begin
gd = 3.0/sqrt(2.)
grad  = [gd 2.0 gd; 2.0 0.0 2.0; gd 2.0 gd]
lap0  = -4.0*2.0- 4.0*gd
#cmax  = 2.0*-3.0/sqrt(2.)
end

# ╔═╡ bd0569a5-0309-4d9d-86ff-e8e6f92d31ec
cs[2, 2]

# ╔═╡ af335a74-58d4-434d-a209-953ecf04aaf5
begin
for i in 1:3
	for j in 1:3
		println("cs [", i,", ", j, "] = ", cs[i, j])
	end
end
end

# ╔═╡ fa190b01-8dec-4d82-9995-96a0afa61be5
begin
coors    = (vec(is), vec(js))
contents = vec(cs)
tsteps    = (1., 1.)
tcl      = jc.clouds(coors, contents, tsteps)
end

# ╔═╡ bf25922a-5b3d-4601-a34b-b014dbd4bd22
reshape(tcl.grad, 3, 3)

# ╔═╡ b6389305-47ff-431d-8990-9fc056c4526e
@assert sum(reshape(tcl.grad, 3, 3) .== reshape(ugrad[aa .> 0], 3, 3)) == 9

# ╔═╡ 227f5f0a-6821-483e-8815-0c4c025f0b31
cs

# ╔═╡ 93f7fe61-de04-4994-8894-4576ed3d701b
begin
xgrad = reshape(tcl.grad, 3, 3) 
end

# ╔═╡ 32fea614-4a42-4cd8-b1ff-58bd53583336
grad

# ╔═╡ 559dd334-de71-4fd1-85c9-43134da59cdf
@assert sum(grad .== xgrad) == 9

# ╔═╡ eb24a7e0-d615-4061-98fa-74ae1ec1ff8e
cs

# ╔═╡ 9d2e8dd5-b4ab-4c51-a2f4-c5d5ad7d077a
xlap = reshape(tcl.lap, 3, 3)

# ╔═╡ 3b8ea9db-b06b-4023-9a7c-911fc4d6eec5
begin
lp11 = (-3.0+3.0)/sqrt(2.)+-2.0+2.0
lp12 = -2.0*2.0/sqrt(2.)-2
lp22 = 4.0*-2.0-4.0*3.0/sqrt(2.)
lap = [lp11 lp12 lp11; lp12 lp22 lp12; lp11 lp12 lp11]
end

# ╔═╡ e3918b67-b2a0-453a-9020-1d9b636d4d0a
@assert sum(lap .== xlap) == 9

# ╔═╡ 05d4215a-3640-40de-8f18-d7109394ed60
begin
cm11 = 2.0/sqrt(2.)
cm12 = 0.
cm22 = -4.
cmax = [cm11 cm12 cm11; cm12 cm22 cm12; cm11 cm12 cm11]
end

# ╔═╡ 8904e52f-6e14-4abb-8042-2307f63746dc
xcmax = reshape(tcl.curmax, 3, 3)

# ╔═╡ 904cf2bb-0dd0-42d1-a418-37d6e549d205
@assert sum(cmax .== xcmax) == 9

# ╔═╡ 3ec495b8-1d51-402f-bbfd-88d2ff7cb23a
begin
cmi11 = -2.0/sqrt(2.)
cmi12 = 0.
cmi22 = -4.0/sqrt(2.)
cmin = [cmi11 cmi12 cmi11; cmi12 cmi22 cmi12; cmi11 cmi12 cmi11]
end

# ╔═╡ 593daba5-16f1-4d73-9228-1198dab656a2
xcmin = reshape(tcl.curmin, 3, 3)

# ╔═╡ f8d614d1-bac5-4d30-a3d1-7a13b8c17e92
cs[1][2]

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

blabel = @bind label Select([:contents, :grad, :lap, :curmax, :curmin, :nodes, :nborders])
	
#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))
	
md"""

Select label to plot $(blabel)

"""
end

# ╔═╡ 7b7981ca-1540-48a1-88e1-4f27e7787b70
md"""
## Graph
"""

# ╔═╡ a779ac6e-5bac-46f1-b8ef-1e3d5b111f4d
md"""

## Code

"""

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
	az, bz, cz = 5., -5., 0.
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

# ╔═╡ 8dca9736-1140-495c-98a3-4cb5acc8ffc1
begin
vals = getfield(xcl, label)
minv, maxv = minimum(vals), maximum(vals)
end;

# ╔═╡ f17d0274-4a61-423c-a76f-870dcef41a60
begin
brange0 = @bind v0 Slider(minv:maxv, default = minv)
brange1 = @bind v1 Slider(minv:maxv, default = maxv)
md"""

Selec range for variable $(label):

minimum $(brange0)
maximum  $(brange1)

"""
end

# ╔═╡ e7544908-23e0-4e3a-ad93-2af5e0dc11f1
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
# ╠═5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╠═a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╠═cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─3922eba2-f322-4b06-b9e0-83bc723d7930
# ╟─3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╠═7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╟─e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╠═9e00d479-f411-4dfd-b51f-95dc9efeecd1
# ╠═a6939579-e310-4ed9-a76d-85d09cba4ce1
# ╠═df582445-6873-4829-9cce-82c1fb737b4b
# ╠═c6a33da6-446f-49bb-acd4-6b024538429e
# ╠═988498e0-6bc9-4c4d-9467-a941a848a157
# ╠═938208b9-9a88-455a-a093-0e9ff4900fe9
# ╠═2c0570a4-c5d0-4a96-9fd9-9618cc1da7f3
# ╠═fcb0a412-ccd6-4406-98f5-cd60ae9f7845
# ╠═2f331294-fa21-431a-928c-0b4d302726ec
# ╠═c70469df-5f59-4ea2-94e0-a247966f0f13
# ╠═792af853-8220-47e8-8a3f-865decbbeafe
# ╠═adf969fd-1555-4099-a197-50281e66597e
# ╠═b8476bc1-8ca1-4539-bdbb-71ddf3fa378c
# ╠═4efe53aa-83b8-46fa-8418-20023d3df966
# ╠═fc36f3ff-0f90-428b-887c-4761706a3dca
# ╠═b4971b5e-a614-4498-a177-ea218ae10933
# ╠═683618f0-6df4-414c-a1d4-6bc27911cc34
# ╠═eb1c838c-42f7-47d9-98d5-9028f834d7c1
# ╠═3638a83e-f986-4752-96e4-a85725387e56
# ╠═7349abcb-1462-43e2-98f5-e4e84d114748
# ╠═17bdc350-b372-40c7-a4b1-6bf5720e03db
# ╠═fc10da1f-8e07-486b-adc5-db1dc4304448
# ╠═d84f1dcc-ba31-4a03-ba60-dcaea8b24938
# ╠═f7035d4a-01aa-4e58-9761-e3e2e346ccd0
# ╠═c66c5247-09de-426b-a9d3-ffb2f69c5edf
# ╠═6b14de6e-b40e-4017-94dd-2d34cfa66000
# ╠═d8626b01-5007-45f2-b47e-b955fc8858ee
# ╠═abda946a-f1b3-4306-9f30-d710fe5e80bd
# ╠═696e4101-c218-4988-b2f5-0d617f386283
# ╠═04b6be28-ed7a-49c6-82c0-1666645e3f98
# ╠═3f03a5ef-ecf4-48dc-bd56-71d365f5e553
# ╠═bbc9921b-8429-40f8-bc02-93b294f3f3d3
# ╠═a1300bcc-ff12-4614-91cf-efbc4791f71a
# ╠═3862afd3-fb51-482e-bbe6-ce8d20ec895b
# ╠═ef73a5ad-e7dc-4bfc-8c06-5680b3f28c85
# ╠═b22dec95-bc58-44fe-b29e-af335c54e710
# ╠═6f5acf59-2a84-4b59-8010-fbd07b769429
# ╠═a0c6aa94-ab32-4d92-a4f9-e17bc7682297
# ╠═f0fe6965-fdbe-4f94-aac3-a41f3cdd63f0
# ╠═d0f3e7c4-d72f-4822-b683-60ba5767278b
# ╠═76f36865-51d6-4eb5-a242-677aa57427e0
# ╠═62e67f7b-bc7a-4b2c-9800-7f89e7ef8d81
# ╠═bf25922a-5b3d-4601-a34b-b014dbd4bd22
# ╠═b01800e5-18b8-4e3a-b1c8-34ada12d9854
# ╠═b6389305-47ff-431d-8990-9fc056c4526e
# ╠═36d0c79d-94b1-4454-b605-50791e51d184
# ╠═13bc2859-55b0-4670-9e42-d42e65d989f0
# ╠═3526ed61-2dc9-49ca-83eb-4feafd29e0ae
# ╠═a796bc34-b937-4f71-a12b-baa5cbd266b9
# ╠═719e52ae-68fe-43c9-877e-1d0e9043e4b7
# ╠═3bbc7b6f-bfa2-4479-bfbc-93dbd695bdbb
# ╠═34923e48-1392-4c29-8c77-f1c58c89a3e2
# ╠═bd0569a5-0309-4d9d-86ff-e8e6f92d31ec
# ╠═af335a74-58d4-434d-a209-953ecf04aaf5
# ╠═fa190b01-8dec-4d82-9995-96a0afa61be5
# ╠═227f5f0a-6821-483e-8815-0c4c025f0b31
# ╠═93f7fe61-de04-4994-8894-4576ed3d701b
# ╠═32fea614-4a42-4cd8-b1ff-58bd53583336
# ╠═559dd334-de71-4fd1-85c9-43134da59cdf
# ╠═eb24a7e0-d615-4061-98fa-74ae1ec1ff8e
# ╠═9d2e8dd5-b4ab-4c51-a2f4-c5d5ad7d077a
# ╠═3b8ea9db-b06b-4023-9a7c-911fc4d6eec5
# ╠═e3918b67-b2a0-453a-9020-1d9b636d4d0a
# ╠═05d4215a-3640-40de-8f18-d7109394ed60
# ╠═8904e52f-6e14-4abb-8042-2307f63746dc
# ╠═904cf2bb-0dd0-42d1-a418-37d6e549d205
# ╠═3ec495b8-1d51-402f-bbfd-88d2ff7cb23a
# ╠═593daba5-16f1-4d73-9228-1198dab656a2
# ╠═f8d614d1-bac5-4d30-a3d1-7a13b8c17e92
# ╠═6c8bf138-8fec-4c69-b4dd-4284faddeed0
# ╟─1a8e9aa9-a47d-40fd-84c6-cfa49f9b1cc4
# ╠═d26c89ae-1629-4e98-8bde-3e8abe8bfd8d
# ╟─5a1832c1-33ff-45dc-8f47-212179dbe862
# ╠═f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
# ╟─4e43c8e3-89e2-44ca-a6ed-48a364d90486
# ╟─13ac9fdf-46d0-4940-80e3-8619f0609108
# ╟─a689debb-8763-45c4-a03d-94c8e970b243
# ╟─8dca9736-1140-495c-98a3-4cb5acc8ffc1
# ╟─f17d0274-4a61-423c-a76f-870dcef41a60
# ╟─e7544908-23e0-4e3a-ad93-2af5e0dc11f1
# ╠═1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
# ╟─7b7981ca-1540-48a1-88e1-4f27e7787b70
# ╟─1c402508-afd3-46a1-8dbc-a23fd9bd63e1
# ╟─a779ac6e-5bac-46f1-b8ef-1e3d5b111f4d
# ╠═d8a02e0a-db35-4965-8322-8741c3ffbd49
# ╠═dfa64554-5fb1-4d63-80d3-19aee7a476b8
# ╠═16b988b0-887f-4672-b347-9c374fcc3fae
