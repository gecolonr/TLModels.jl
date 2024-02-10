using Parameters

@with_kw mutable struct LineModelParams
    z_km::Union{Vector{ComplexF64}, ComplexF64}
    y_km::Union{Vector{ComplexF64}, ComplexF64}
    z_km_ω::ComplexF64
    Z_c::Float64
    M::Int64
    l_dict::Union{Dict{String, Int64}, Nothing} = nothing
    alg_line_name::String
    l_seg::Union{Int64, Float64}
    line_scale::Float64
    load_scale::Float64
end

default_2_bus_line_dict = Dict(
    "BUS 1-BUS 2-i_1" => 100,
    "BUS 1-BUS 2-i_1_static" => 100,
)

default_9_bus_line_dict = Dict(
    "Bus 5-Bus 4-i_1" => 90,
    "Bus 7-Bus 8-i_1" => 80,
    "Bus 6-Bus 4-i_1" => 100,
    "Bus 7-Bus 5-i_1" => 170,
    "Bus 8-Bus 9-i_1" => 110,
    "Bus 9-Bus 6-i_1" => 180,
)

export LineModelParams

export default_2_bus_line_dict
export default_9_bus_line_dict