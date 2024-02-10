module TLmodels

using CSV
using DataFrames
using PowerSystems

include("TLmodels_structs.jl")

function get_line_parameters_from_data(impedance_csv, capacitance_csv, M)
    # Data needs to be imported in a particular format. Each row is the value of M, and 
    # each column is the value of R and L. For more M >= 1, the row will have more columns
    # corresponding to the branches.
    
    df_imp = CSV.read(impedance_csv, DataFrame);
    c_km = CSV.read(capacitance_csv, DataFrame)."C"

    r_km = vec(zeros(M, 1))
    l_km = vec(zeros(M, 1))
    
    for m in 1 : M
        r_km[m,1] = df_imp[M, "R"*string(m)]
        l_km[m,1] = df_imp[M, "L"*string(m)]
    end
    
    # We nominally use 60 Hz as the grid frequency.
    f = 60
    ω = 2*pi*f

    # Compute impedance of each branch of the line and store in the M dimensional vector z_km
    # Compute the admittance of each line and store in the scalar y_km
    x_km = ω*l_km
    z_km = r_km + im*x_km
    g_km = [0.0]
    b_km = ω*c_km
    y_km = (g_km + im*b_km)[1]
    
    # Compute the equivalent impedance of the line at 60 Hz and its characteristic impedance
    # for normalizing the line parameters
    Y_ = 0
    for i in 1:M
        Y_ = Y_ + 1/(z_km[i])
    end
    z_km_ω = 1/Y_
    Z_c = abs(sqrt(z_km_ω/y_km))
    
    return z_km, y_km, z_km_ω, Z_c 
end

function create_statpi_system(sys, p::LineModelParams)

    sys_statpi = deepcopy(sys)

    # Unpacking parameters 
    z_km_ω = p.z_km_ω # Ω/km
    y_km = p.y_km # S/km
    Z_c = p.Z_c # Ω
    
    z_km_ω_pu = z_km_ω/Z_c
    y_km_pu = y_km*Z_c
    γ = sqrt(z_km_ω*y_km)
    
    line_scale = p.line_scale

    # Rewrite all default line impedances with impednaces from the data
    for ll in get_components(Line, sys_statpi)
        l = p.l_dict[ll.name]*line_scale #km
        z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
        y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
        ll.r = real(z_ll)
        ll.x = imag(z_ll)
        ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
    end

    return sys_statpi
end

function create_dynpi_system(sys, p::LineModelParams)

    sys_dynpi = deepcopy(sys)

    # Unpacking parameters 
    z_km_ω = p.z_km_ω # Ω/km
    y_km = p.y_km # S/km
    Z_c = p.Z_c # Ω
    
    z_km_ω_pu = z_km_ω/Z_c
    y_km_pu = y_km*Z_c
    γ = sqrt(z_km_ω*y_km)
    
    line_scale = p.line_scale

    # Rewrite all default line impedances with impednaces from the data
    for ll in get_components(Line, sys_dynpi)
        l = p.l_dict[ll.name]*line_scale #km
        z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
        y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
        ll.r = real(z_ll)
        ll.x = imag(z_ll)
        ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
        # Dynamicize all lines except for one and add to the system
        if (ll.name != p.alg_line_name)
            dyn_branch = DynamicBranch(get_component(Line, sys_dynpi, ll.name))
            add_component!(sys_dynpi, dyn_branch)
        end
    end

    return sys_dynpi
end

function create_MSSB_system(sys, p::LineModelParams)
    
    sys_MSSB = deepcopy(sys)
    # Unpacking parameters 
    z_km_ω = p.z_km_ω # Ω/km
    y_km = p.y_km # S/km
    Z_c = p.Z_c # Ω
    
    z_km_ω_pu = z_km_ω/Z_c
    y_km_pu = y_km*Z_c
    γ = sqrt(z_km_ω*y_km)
    
    line_scale = p.line_scale

    # How many line segments to divide each line into.
    l_seg = p.l_seg
    
    for ll in collect(get_components(Line, sys_MSSB))
        
        # Do not segment or dynamicize the line to be tripped
        if ll.name == p.alg_line_name
            l = p.l_dict[ll.name]*line_scale #km
            z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
            y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
            ll.r = real(z_ll)
            ll.x = imag(z_ll)
            ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
            continue
        end
        
        l = p.l_dict[ll.name]*line_scale #km
        N = Int(ceil(l/l_seg))
        # l_seg = l/N

        z_km_ω_pu_seg = z_km_ω_pu*l_seg
        y_km_pu_seg = y_km_pu*l_seg

        bus_from = ll.arc.from
        bus_to = ll.arc.to
        
        # Create a bunch of Bus
        start_bus = bus_from
        for b_ix in 1 : N - 1
            # In a previous version of PSY, ACBus, was simply called Bus. Verify you are using latest PSY version.
            bus_to_create = ACBus(
                number = 1000000000 + 100000*bus_from.number + 100*bus_to.number + b_ix,
                name = bus_from.name * "-" * bus_to.name * "-internal-bus_" * string(b_ix),
                bustype = ACBusTypes.PQ,
                angle = bus_from.angle,
                magnitude = bus_from.magnitude,
                voltage_limits = bus_from.voltage_limits,
                base_voltage = bus_from.base_voltage,
                area = bus_from.area,
                load_zone = bus_from.load_zone
            )
            add_component!(sys_MSSB, bus_to_create)
            end_bus = bus_to_create

            line_to_create = Line(
                name = ll.name * "_segment_" * string(b_ix),
                available = true,
                active_power_flow = ll.active_power_flow,
                reactive_power_flow = ll.reactive_power_flow,
                arc = Arc(from = start_bus, to = end_bus),
                r = real(z_km_ω_pu_seg),
                x = imag(z_km_ω_pu_seg),
                b = (from = imag(y_km_pu_seg)/(2), to = imag(y_km_pu_seg)/(2)),
                rate = ll.rate,
                angle_limits = ll.angle_limits,
            )
            add_component!(sys_MSSB, line_to_create)
            dyn_branch = DynamicBranch(get_component(Line, sys_MSSB, line_to_create.name))
            add_component!(sys_MSSB, dyn_branch)
            
            start_bus = end_bus
        end

        line_to_create = Line(
                name = ll.name * "_segment_" * string(N),
                available = true,
                active_power_flow = ll.active_power_flow,
                reactive_power_flow = ll.reactive_power_flow,
                arc = Arc(from = start_bus, to = bus_to),
                r = real(z_km_ω_pu_seg),
                x = imag(z_km_ω_pu_seg),
                b = (from = imag(y_km_pu_seg)/(2), to = imag(y_km_pu_seg)/(2)),
                rate = ll.rate,
                angle_limits = ll.angle_limits,
        )

        add_component!(sys_MSSB, line_to_create)
        dyn_branch = DynamicBranch(get_component(Line, sys_MSSB, line_to_create.name))
        add_component!(sys_MSSB, dyn_branch)

        remove_component!(sys_MSSB, ll)
    end
    
    return sys_MSSB
end

function create_MSMB_system(sys, p::LineModelParams)
    
    sys_MSMB = deepcopy(sys)
   
    # Unpacking parameters 
    z_km = p.z_km # Ω/km
    z_km_ω = p.z_km_ω # Ω/km
    y_km = p.y_km # S/km
    Z_c = p.Z_c # Ω
    
    z_km_pu = z_km/Z_c
    z_km_ω_pu = z_km_ω/Z_c
    y_km_pu = y_km*Z_c
    γ = sqrt(z_km_ω*y_km)
    
    line_scale = p.line_scale

    # How many line segments to divide each line into.
    M = p.M
    l_seg = p.l_seg
    
    for ll in collect(get_components(Line, sys_MSMB))
        
        # Do not segment or dynamicize the line to be tripped
        if ll.name == p.alg_line_name
            l = p.l_dict[ll.name]*line_scale #km
            z_ll = z_km_ω_pu*l*(sinh(γ*l)/(γ*l))
            y_ll = y_km_pu*l*(tanh(γ*l/2)/(γ*l/2))
            ll.r = real(z_ll)
            ll.x = imag(z_ll)
            ll.b = (from = imag(y_ll)/2, to = imag(y_ll)/2)
            continue
        end
        
        l = p.l_dict[ll.name]*line_scale #km
        N = Int(ceil(l/l_seg))
        # l_seg = l/N

        z_km_pu_seg = z_km_pu*l_seg
        y_km_pu_seg = y_km_pu*l_seg

        bus_from = ll.arc.from
        bus_to = ll.arc.to
        
        # Create a bunch of Bus
        start_bus = bus_from
        for b_ix in 1 : N - 1

            bus_to_create = ACBus(
                number = 1000000000 + 100000*bus_from.number + 100*bus_to.number + b_ix,
                name = bus_from.name * "-" * bus_to.name * "-internal-bus_" * string(b_ix),
                bustype = ACBusTypes.PQ,
                angle = bus_from.angle,
                magnitude = bus_from.magnitude,
                voltage_limits = bus_from.voltage_limits,
                base_voltage = bus_from.base_voltage,
                area = bus_from.area,
                load_zone = bus_from.load_zone
            )

            add_component!(sys_MSMB, bus_to_create)
            end_bus = bus_to_create
           
            for m in 1 : M
                line_to_create = Line(
                    name = ll.name * "_segment_" * string(b_ix) * "_branch_" * string(m),
                    available = true,
                    active_power_flow = ll.active_power_flow,
                    reactive_power_flow = ll.reactive_power_flow,
                    arc = Arc(from = start_bus, to = end_bus),
                    r = real(z_km_pu_seg[m]),
                    x = imag(z_km_pu_seg[m]),
                    b = (from = imag(y_km_pu_seg)/(2*M), to = imag(y_km_pu_seg)/(2*M)),
                    rate = ll.rate,
                    angle_limits = ll.angle_limits,
                )
                add_component!(sys_MSMB, line_to_create)
                dyn_branch = DynamicBranch(get_component(Line, sys_MSMB, line_to_create.name))
                add_component!(sys_MSMB, dyn_branch)
            end

            start_bus = end_bus
        end
        
        for m in 1 : M
            line_to_create = Line(
                    name = ll.name * "_segment_" * string(N) * "_branch_" * string(m),
                    available = true,
                    active_power_flow = ll.active_power_flow,
                    reactive_power_flow = ll.reactive_power_flow,
                    arc = Arc(from = start_bus, to = bus_to),
                    r = real(z_km_pu_seg[m]),
                    x = imag(z_km_pu_seg[m]),
                    b = (from = imag(y_km_pu_seg)/(2*M), to = imag(y_km_pu_seg)/(2*M)),
                    rate = ll.rate,
                    angle_limits = ll.angle_limits,
            )

            add_component!(sys_MSMB, line_to_create)
            dyn_branch = DynamicBranch(get_component(Line, sys_MSMB, line_to_create.name))
        add_component!(sys_MSMB, dyn_branch)
        end
        remove_component!(sys_MSMB, ll)
    end
    
    return sys_MSMB
end

function create_all_TL_model_systems(sys, line_params)

    sys_statpi = create_statpi_system(sys, line_params) 
    sys_dynpi = create_dynpi_system(sys, line_params)
    sys_MSSB = create_MSSB_system(sys, line_params)
    sys_MSMB = create_MSMB_system(sys, line_params)

    return sys_statpi, sys_dynpi, sys_MSSB, sys_MSMB
end

export get_line_parameters_from_data
export create_statpi_system, create_dynpi_system 
export create_MSSB_system, create_MSMB_system
export create_all_TL_model_systems

end
