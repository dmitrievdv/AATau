module AATau

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

end

function screenτ(screen::DustScreen, Δh::Real)
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

function BV(star::AATauStar, i::Real, screen::DustScreen, yₛ::Real; h = 0.02, scatter = 1, sym = true)
  coords = [-1:h:1;]
  EV = 0
  EB = 0
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

    # print(τ, " ", y - yₛ, "\n")
    cont = spotorstar(star, i, x, y)
    ΔV = τ/log(2.512)
    ΔB = ΔV/screen.Rᵥ + ΔV
    EV += cont(551e-7)*90e-7*(exp(-τ) + 0.125*scatter)
    EB += cont(445e-7)*94e-7*(2.512^(-ΔB) + 0.175*scatter)
  end
  EV *= h^2*star.R^2
  EB *= h^2*star.R^2
  return -2.5*log10(EB), -2.5*log10(EV)
end

function BV(star::AATauStar, i::Real, screen::DustScreen, yₛ::Array{T, 1}; h = 0.02, scatter = 1, sym = true) where {T<:Real}
  BVss(i::Real, y::Real) = BV(star, i, screen, y, h = h, scatter = scatter, sym = sym)
  B0, V0 = BVss(i, -1e99)
  return first.(BVss.(i, yₛ)) .- B0, last.(BVss.(i, yₛ)) .- V0
end

end