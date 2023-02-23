using Plots

mutable struct Universe
	Ω_m::Float64
	Ω_k::Float64
	Ω_Λ::Float64
end

einsteinDeSitter = Universe(1.0, 0.0, 0.0)
flatUniverse1 = Universe(0.3, 0.0, 0.7)

H_0 = 70.0	# km/s/Mpc
c = 1.0 # 299792.458	# km/s
L = 1.0		# Mpc

zpoints = (0.1:0.01:10)

function f(z::Float64)
	return 1 / sqrt(universe.Ω_m * (1+z)^3 
		+ universe.Ω_k * (1+z)^2 + universe.Ω_Λ)
end

function integrate_simpson(a::Float64, b::Float64)
	N = 10000
	h = (b-a)/N

	I = f(a) + f(b)
	for i in 1:(N-1)
		I += 4 * f(a + i*h)
	end
	for i in 2:(N-2)
		I += 2 * f(a + i*h)
	end
	I *= (1/3) * h
	return I
end


function get_angular_scale(universe)
	universe = universe

	θpoints = zeros(length(zpoints))
	dApoints = zeros(length(zpoints))

	for i in 1:length(zpoints)
		z = zpoints[i]
		χ = (c/H_0) * integrate_simpson(0.0,z)
		a = 1/(1+z)
		dApoints[i] = a * χ
		θpoints[i] = L / (a * χ)
	end

	θpoints = θpoints ./ L

	return [θpoints, dApoints]
end

function r()
	include("angularsize.jl")
end

universe = einsteinDeSitter
edsdata = get_angular_scale(einsteinDeSitter)
universe = flatUniverse1
fu1data = get_angular_scale(flatUniverse1)

println("Ω_m = ", universe.Ω_m, "\nΩ_k = ", universe.Ω_k, "\nΩ_Λ = ", universe.Ω_Λ)

function plotθ()
	plot(size=(600,400), titlefontsize=28, font=("monospace", 28))
	plot!(zpoints, edsdata[1], 
		label="einstein-de sitter: Ω_m = 1.0, Ω_k = 0.0, Ω_Λ = 0.0")
	plot!(zpoints, fu1data[1], 
		label="flat universe 1:    Ω_m = 0.3, Ω_k = 0.0, Ω_Λ = 0.7")
	xticks!(0:1:10)
	xlabel!("redshift, z")
	ylabel!("θ(z)/L")
	png("plot_θ_z.png")
end

function plotdA()
	plot(size=(600,400), titlefontsize=28, font=("monospace", 28))
	plot!(zpoints, edsdata[2], label="einstein-de sitter, Ω_m = 1.0, Ω_k = 0.0, Ω_Λ = 0.0")
	plot!(zpoints, fu1data[2], label="flat universe 1, Ω_m = 0.3, Ω_k = 0.0, Ω_Λ = 0.7")
	xticks!(0:1:10)
	xlabel!("redshift, z")
	ylabel!("d_A (Mpc)")
	png("plot_dA_z.png")
end

plotθ()
plotdA()
