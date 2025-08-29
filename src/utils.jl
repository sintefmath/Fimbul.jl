function thermal_radius_aquifer(Vin, Haq, phi, Cf, Cr)

    ϕ = phi
    Caq = Cf*ϕ + Cr*(1 - ϕ)
    return sqrt(Cf*Vin/(Caq*π*Haq))

end