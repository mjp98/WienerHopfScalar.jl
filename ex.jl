### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ bd829ece-2476-11ed-1d11-8bc94bcf2111
begin 
	using Pkg
	Pkg.activate(".scalarwienerhopf")
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

# ╔═╡ 2c3e451d-23d6-43f9-931d-627ccb807576
Pkg.add(PackageSpec(url="https://github.com/mjp98/Aeroacoustics.jl"));

# ╔═╡ 80d3ed34-12e0-4c13-ba04-34c2b875a7f0
begin 
	using WienerHopfScalar
	using HolomorphicFun
	using PhaseColors
	using Plots
end

# ╔═╡ d9fce9ca-d87f-4e2a-99ae-bf7b74a80fcc
using Aeroacoustics

# ╔═╡ 996d9507-b396-4b08-ba81-5287dc840c00
md"""

# Scalar Wiener Hopf

WienerHopfScalar.jl

"""

# ╔═╡ 8112ea74-7aed-40db-a944-a1453fab43ed
HolomorphicFun.ScalarConstant(0.2)

# ╔═╡ ab3c9340-350e-4baf-838b-5b34679ff6a2
begin
	let 
		z = PhaseColors.ℂ(5;n=500);
		x,y = PhaseColors.ℂ2xy(z);
		K = WienerHopfScalar.GammaKernel(0.4);

		plot(x,y, PhaseColors.portrait(K.(z),PhaseColors.blueorange),yflip=false)
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
	
	evaluate(K::TestKernel, z) = 2 + (z+1)/sqrt(im*(z-1))/sqrt(-im*(z+1))
	
	leftlimit(K::TestKernel) = 1
	rightlimit(K::TestKernel) = 3
	
	@nopoles TestKernel
	@noroots TestKernel
end

# ╔═╡ 7427ac26-5586-4b53-8235-278ae7ededab
f(z) = sqrt(im*(z-1))*sqrt(-im*(z+1))

# ╔═╡ e4a588fa-edb7-4c5e-8b68-2d45ae1aa474
f(1e6)

# ╔═╡ 295ae430-6174-4a4a-b9fc-c236a03648ac
begin
	let
		z = PhaseColors.ℂ(5);
		x,y = PhaseColors.ℂ2xy(z);
		K = TestKernel(0.4);

		plot(x,y, PhaseColors.portrait(K.(z),PhaseColors.math),yflip=false)
	end
end

# ╔═╡ 098700b2-2bf0-49ca-a21f-941cf3246aee
begin
	
	let
		K = TestKernel(0.4);
		L = factorise(K)

		z = PhaseColors.ℂ(5);
		x,y = PhaseColors.ℂ2xy(z);

		plot(x,y, PhaseColors.portrait(L.(z,false),PhaseColors.:math),yflip=false)
		
	end
end

# ╔═╡ 0958fc8a-d8ed-4119-886c-6fa0185a5a15
spectrum = GaussianEnergySpectrum(0.2)

# ╔═╡ Cell order:
# ╟─996d9507-b396-4b08-ba81-5287dc840c00
# ╠═bd829ece-2476-11ed-1d11-8bc94bcf2111
# ╠═791f2892-0df1-4127-9b10-3f545c100e0c
# ╠═80d3ed34-12e0-4c13-ba04-34c2b875a7f0
# ╠═8112ea74-7aed-40db-a944-a1453fab43ed
# ╠═ab3c9340-350e-4baf-838b-5b34679ff6a2
# ╠═552e513a-75b9-4a94-b7c0-19b876fa0b9c
# ╠═baeaaf12-10e1-49f6-a890-d6e62aba1b8e
# ╠═7427ac26-5586-4b53-8235-278ae7ededab
# ╠═e4a588fa-edb7-4c5e-8b68-2d45ae1aa474
# ╠═295ae430-6174-4a4a-b9fc-c236a03648ac
# ╠═098700b2-2bf0-49ca-a21f-941cf3246aee
# ╠═2c3e451d-23d6-43f9-931d-627ccb807576
# ╠═d9fce9ca-d87f-4e2a-99ae-bf7b74a80fcc
# ╠═0958fc8a-d8ed-4119-886c-6fa0185a5a15
