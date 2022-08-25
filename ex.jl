### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ bd829ece-2476-11ed-1d11-8bc94bcf2111
begin 
	using Pkg
	Pkg.activate(".new")
		Pkg.add(PackageSpec(url="https://github.com/mjp98/HolomorphicFun.jl"));
		Pkg.add(PackageSpec(url="https://github.com/mjp98/RiemannHilbert.jl",rev="update1_8"));
		Pkg.add(PackageSpec(url="https://github.com/mjp98/WienerHopf.jl",rev="update1_8"));
		Pkg.add(PackageSpec(url="https://github.com/mjp98/ScaledLines.jl"));
		Pkg.add(PackageSpec(url="https://github.com/mjp98/WienerHopfScalar.jl"));
end;

# ╔═╡ 791f2892-0df1-4127-9b10-3f545c100e0c
begin
	Pkg.add(PackageSpec(url="https://github.com/mjp98/PhaseColors.jl"));
	Pkg.add("Plots");
end;

# ╔═╡ 80d3ed34-12e0-4c13-ba04-34c2b875a7f0
begin 
	using WienerHopfScalar
	using HolomorphicFun
	using PhaseColors
	using Plots
end

# ╔═╡ 8112ea74-7aed-40db-a944-a1453fab43ed
HolomorphicFun.ScalarConstant(0.2)

# ╔═╡ ab3c9340-350e-4baf-838b-5b34679ff6a2
begin
	let 
		z = PhaseColors.ℂ(5);
		x,y = PhaseColors.ℂ2xy(z);
		K = WienerHopfScalar.GammaKernel(0.4);

		plot(x,y, PhaseColors.portrait(K.(z),PhaseColors.math),yflip=false)
	end
end

# ╔═╡ 552e513a-75b9-4a94-b7c0-19b876fa0b9c
md"## Define a new kernel"

# ╔═╡ baeaaf12-10e1-49f6-a890-d6e62aba1b8e
begin
	import WienerHopfScalar: evaluate, leftlimit, rightlimit
	import WienerHopfScalar: candidate_poles, candidate_roots
	
	struct TestKernel{T} <: WienerHopfKernel
	    μ::T
	end
	
	evaluate(K::TestKernel, z) = 1 + K.μ/sqrt(im*(z-1))/sqrt(-im*(z+1))
	
	leftlimit(K::TestKernel) = 1
	rightlimit(K::TestKernel) = 1
	
	@nopoles TestKernel
	@noroots TestKernel
end

# ╔═╡ 7427ac26-5586-4b53-8235-278ae7ededab


# ╔═╡ 295ae430-6174-4a4a-b9fc-c236a03648ac
begin
	let
		z = PhaseColors.ℂ(5);
		x,y = PhaseColors.ℂ2xy(z);
		K = TestKernel(0.4);

		plot(x,y, PhaseColors.portrait(K.(z),PhaseColors.math),yflip=false)
	end
end

# ╔═╡ 58d7aa6b-a7f1-4444-a54b-9445f28793ee
c

# ╔═╡ 098700b2-2bf0-49ca-a21f-941cf3246aee
begin
	
	let
		K = TestKernel(0.4);
		L = factorise(K)

		z = PhaseColors.ℂ(5);
		x,y = PhaseColors.:ℂ2xy(z);

		plot(x,y, PhaseColors.portrait(L.(z,true),PhaseColors.:math),yflip=false)
		
	end
end

# ╔═╡ Cell order:
# ╠═bd829ece-2476-11ed-1d11-8bc94bcf2111
# ╠═791f2892-0df1-4127-9b10-3f545c100e0c
# ╠═80d3ed34-12e0-4c13-ba04-34c2b875a7f0
# ╠═8112ea74-7aed-40db-a944-a1453fab43ed
# ╠═ab3c9340-350e-4baf-838b-5b34679ff6a2
# ╠═552e513a-75b9-4a94-b7c0-19b876fa0b9c
# ╠═baeaaf12-10e1-49f6-a890-d6e62aba1b8e
# ╠═7427ac26-5586-4b53-8235-278ae7ededab
# ╠═295ae430-6174-4a4a-b9fc-c236a03648ac
# ╠═58d7aa6b-a7f1-4444-a54b-9445f28793ee
# ╠═098700b2-2bf0-49ca-a21f-941cf3246aee
