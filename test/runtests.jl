cd(@__DIR__)
using TLmodels
using CSV
using DataFrames
using PowerSystems
using Revise
using Test


const PSY = PowerSystems

impedance_csv = "../../../../TLModels/data/cable_data/dommel_data.csv"
capacitance_csv = "../../../../TLModels/data/cable_data/dommel_data_C.csv"

M = 3
z_km, y_km, z_km_ω, Z_c = get_line_parameters_from_data(impedance_csv, capacitance_csv, M)

line_params = LineModelParams(
    z_km, 
    y_km, 
    z_km_ω, 
    Z_c,
    M,
    default_9_bus_line_dict,    
    "BUS 5-BUS 4-i_1",
    10.0,
    1.0,
    1.0
)  

file_name = "../../../../TLModels/data/json_data/9bus_VSM_SM_GFL_.json"

sys = System(joinpath(pwd(), file_name));

sys_statpi, sys_dynpi, sys_MSSB, sys_MSMB = create_all_TL_model_systems(sys, line_params)

using PowerFlows
sol_statpi = solve_powerflow(ACPowerFlow(), sys_statpi)
sol_statpi["bus_results"]

sol_dynpi = solve_powerflow(ACPowerFlow(), sys_dynpi)
sol_dynpi["bus_results"]

sol_MSSB = solve_powerflow(ACPowerFlow(), sys_MSSB)
sol_MSSB["bus_results"]

sol_MSMB = solve_powerflow(ACPowerFlow(), sys_MSMB)
sol_MSMB["bus_results"]

@testset "TLmodels.jl" begin
    # Write your tests here.
end
