module AATau

using Printf
using Plots
# using PyPlot

const R☉ = 6.955e10
const M☉ = 1.989e33
const year_seconds = 3.155e7
const G = 6.67259e-8
const σ = 5.67e-5


struct AATauStar
  R::Float64
  M::Float64
  T::Float64
  Tspot::Float64
  Mdot::Float64
  rmi::Float64
  rmo::Float64
  vesc::Float64
  starcont
  spotcont
  function AATauStar(R, M, T, Mdot, rmi, rmo)
    R = R*R☉; M = M*M☉
    Mdot = Mdot*M☉/year_seconds
    vesc = √(2G*M/R)
    h1 = asin(√(1/rmo)); h2 = asin(√(1/rmi))
    Tspot = (vesc^2/2*Mdot*(1 - 2/(1/sin(h1)^2+1/sin(h2)^2))/R^2/4/π/abs(cos(h2)-cos(h1))/σ)^0.25
    if Tspot < T
      Tspot = T
    end
    new(R, M, T, Tspot, Mdot, rmi*R, rmo*R, vesc, planck(T), planck(Tspot))
  end
end

struct DustScreen
  h::Real
  τ::Real
  Rᵥ::Real
  Rᵤ::Real
  Rᵢ::Real
  DustScreen(h, τ; Rᵥ = 3.1, Rᵤ = 4.1, Rᵢ = 3.2) = new(h, τ, 3.1, 4.1, 3.2)
end

function screenτ(screen::DustScreen, Δh::Real)
  if(screen.h ≈ 0) 
  	return 0
  end
  return screen.τ*exp(-abs(Δh/screen.h))
end

function planck(T::Real)
  c = 2.99792458e10
  h = 6.626176e-27
  k = 1.380649e-16

  return f(λ::Real) = 2*h/c^2/λ^5 * 1.0/(exp(h*c/λ/(k*T))-1)
end

function spotorstar(star::AATauStar, i::Real,  x::Real, y::Real)
  θ = i/180*π
  inv_rm = sin(θ)^2 + y^2*cos(2*θ) + x^2*cos(θ)^2 - y*√(1-x^2-y^2)*sin(2*θ)
  inv_rmo = star.R/star.rmo; inv_rmi = star.R/star.rmi
  if(inv_rmo <= inv_rm <= inv_rmi)
    # print("hot!\n")
    return star.spotcont
  else
    # print("cold!\n")
    return star.starcont
  end
end

function UBVI(star::AATauStar, i::Real, screen::DustScreen, yₛ::Real; h = 0.02, scatter = 1, sym = true)
  coords = [-1:h:1;]
  EV = 0
  EB = 0
  EU = 0; EI = 0
  for x in coords, y in coords
    # print(x, " ", y, " ", (y < yₛ) & !sym, "\n")
    r² = x^2+y^2
    if r² > 0.999999
      continue
    end
    # print(r², "\n")
    τ = screenτ(screen, y - yₛ)
    if ((y < yₛ) && !sym)
      # print("low\n")
      τ = screen.τ
    end
    
    # println(τ, " ", y - yₛ)
    cont = spotorstar(star, i, x, y)
    ΔV = τ/log(2.512)
    ΔB = ΔV/screen.Rᵥ + ΔV
    ΔU = ΔV/screen.Rᵤ + ΔB
    ΔI = ΔV/screen.Rᵢ + ΔV
    EV += cont(551e-7)*90e-7*(exp(-τ) + 0.121*scatter)
    EB += cont(445e-7)*94e-7*(2.512^(-ΔB) + 0.167*scatter)
    EU += cont(365e-7)*66e-7*(2.512^(-ΔU) + 0.222*scatter)
    EI += cont(806e-7)*149e-7*(2.512^(-ΔI) + 0.065*scatter)
  end
  EV *= h^2*star.R^2
  EB *= h^2*star.R^2
  EU *= h^2*star.R^2
  EI *= h^2*star.R^2
  return -2.5*log10(EU), -2.5*log10(EB), -2.5*log10(EV), -2.5*log10(EI)
end

function UBVI(star::AATauStar, i::Real, screen::DustScreen, yₛ::Array{T, 1}; h = 0.02, scatter = 1, sym = true) where {T<:Real}
  UBVIss(i::Real, y::Real) = UBVI(star, i, screen, y, h = h, scatter = scatter, sym = sym)
  U0, B0, V0, I0 = UBVIss(i, -1e99)
  UBVIt = UBVIss.(i, yₛ)
  U = [UBVIp[1] for UBVIp in UBVIt]
  B = [UBVIp[2] for UBVIp in UBVIt]
  V = [UBVIp[3] for UBVIp in UBVIt]
  I = [UBVIp[4] for UBVIp in UBVIt]
  return Dict("U" => U .- U0, "B" => B .- B0, "V" => V .- V0, "I" => I .- I0)
end

function saveviscontinuummap(star::AATauStar, i::Real, screen::DustScreen, yₛ; h = 0.02, sym = true)
  coords = [-1:h:1;]
  open("vis_continuum_map.out", "w") do io
    for x in coords
      for y in coords
        r² = x^2 + y^2
        if r² > 0.999999
          @printf(io, "%8.f %8.f %13.e\n", x, y, 0.0)
          continue
        end
    
        τ = screenτ(screen, y - yₛ)
        if ((y < yₛ) && !sym)
          # print("low\n")
          τ = screen.τ
        end
        cont = spotorstar(star, i, x, y)
        @printf(io, "%8.f %8.f %13.e\n", x, y, cont(551e-7)*90e-7*exp(-τ))
      end
      print(io, "\n")
    end
  end
end

function manyscreens(star, i, screens, ys; sym = false, scatter = 1, outfile = "manyscreens.dat")
  plt = Plots.scatter([0], [0], yflip = true, ylabel = "V", xlabel = "B-V", label = "Uneclipsed")
  f = open(outfile, "w")
  for screen in screens 
    UBVI = AATau.UBVI(star, i, screen, ys, sym = sym, scatter = scatter)
    plot!(plt, UBVI["B"].-UBVI["V"], UBVI["V"], label = string(screen.h))
    # plot_title!(plt, "̇M = "*string(star.Mdot)*"; mag = "*string(star.rmi)*"-"*string(star.rmo))
    println(f, "#h= ", screen.h)
    println(f, 0.0, " ", 0.0, " ", 0.0, " ", 0.0, " ", screen.h)
    for (u,b,v,i) in zip(UBVI["U"], UBVI["B"], UBVI["V"], UBVI["I"])
    	println(f, u, " ", b, " ", v, " ", i, " ", screen.h)
    end
    println(f, "")
  end
  close(f)
  return plt
end
end