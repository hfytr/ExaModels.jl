using ExaModels, ExaPowerIO

function parse_ac_power_struct_data(filename, backend)
    data = ExaPowerIO.parse_pglib(Float32, filename, "../data"; out_type=ExaPowerIO.Data)
    convert = a -> ExaModels.convert_array(a, backend)
    arc_from = [(i, b.fbus, b.tbus) for (i, b) in enumerate(data.branch)]
    arc_to = [(i, b.tbus, b.fbus) for (i, b) in enumerate(data.branch)]
    arc = [arc_from; arc_to]
    (
        bus = convert(data.bus),
        gen = convert(data.gen),
        arc = convert(arc),
        arc_idxs = convert([(fidx=i, tidx=i+length(data.branch)) for i in 1:length(data.branch)]),
        branch = convert(data.branch),
        ref_buses = convert([i for i in 1:length(data.bus) if data.bus[i].type == 3]),
        vmax = convert([bu.vmax for bu in data.bus]),
        vmin = convert([bu.vmin for bu in data.bus]),
        pmax = convert([g.pmax for g in data.gen]),
        pmin = convert([g.pmin for g in data.gen]),
        qmax = convert([g.qmax for g in data.gen]),
        qmin = convert([g.qmin for g in data.gen]),
        angmax = convert([br.angmax for br in data.branch]),
        angmin = convert([br.angmin for br in data.branch]),
        rate_a = convert([data.branch[l].ratea for (l, i, j) in arc]),
    )
end

function _exa_ac_power_struct_model(backend, filename)

    data = parse_ac_power_struct_data(filename, backend)

    w = ExaModels.ExaCore(backend = backend)

    va = ExaModels.variable(w, length(data.bus);)

    vm = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = data.vmin,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        w,
        g.c[1] * pg[i]^2 + g.c[2] * pg[i] + g.c[3] for (i, g) in enumerate(data.gen)
    )

    c1 = ExaModels.constraint(w, va[b] for b in data.ref_buses)

    c2 = ExaModels.constraint(
        w,
        p[i.fidx] - b.c5 * vm[b.fbus]^2 -
        b.c3 * (vm[b.fbus] * vm[b.tbus] * cos(va[b.fbus] - va[b.tbus])) -
        b.c4 * (vm[b.fbus] * vm[b.tbus] * sin(va[b.fbus] - va[b.tbus]))
        for (i, b) in zip(data.arc_idxs, data.branch)
    )

    c3 = ExaModels.constraint(
        w,
        q[i.fidx] +
        b.c6 * vm[b.fbus]^2 +
        b.c4 * (vm[b.fbus] * vm[b.tbus] * cos(va[b.fbus] - va[b.tbus])) -
        b.c3 * (vm[b.fbus] * vm[b.tbus] * sin(va[b.fbus] - va[b.tbus]))
        for (i, b) in zip(data.arc_idxs, data.branch)
    )

    c4 = ExaModels.constraint(
        w,
        p[i.tidx] - b.c7 * vm[b.tbus]^2 -
        b.c1 * (vm[b.tbus] * vm[b.fbus] * cos(va[b.tbus] - va[b.fbus])) -
        b.c2 * (vm[b.tbus] * vm[b.fbus] * sin(va[b.tbus] - va[b.fbus]))
        for (i, b) in zip(data.arc_idxs, data.branch)
    )

    c5 = ExaModels.constraint(
        w,
        q[i.tidx] +
        b.c8 * vm[b.tbus]^2 +
        b.c2 * (vm[b.tbus] * vm[b.fbus] * cos(va[b.tbus] - va[b.fbus])) -
        b.c1 * (vm[b.tbus] * vm[b.fbus] * sin(va[b.tbus] - va[b.fbus]))
        for (i, b) in zip(data.arc_idxs, data.branch)
    )

    c6 = ExaModels.constraint(
        w,
        va[b.fbus] - va[b.tbus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )
    c7 = ExaModels.constraint(
        w,
        p[i.fidx]^2 + q[i.fidx]^2 - b.ratea^2 for (i, b) in zip(data.arc_idxs, data.branch);
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        w,
        p[i.tidx]^2 + q[i.tidx]^2 - b.ratea^2 for (i, b) in zip(data.arc_idxs, data.branch);
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * vm[i]^2 for (i, b) in enumerate(data.bus))
    c10 = ExaModels.constraint(w, b.qd - b.bs * vm[i]^2 for (i, b) in enumerate(data.bus))

    c11 = ExaModels.constraint!(w, c9, a_bus => p[i] for (i, (l, a_bus, j)) in enumerate(data.arc))
    c12 = ExaModels.constraint!(w, c10, a_bus => q[i] for (i, (l, a_bus, j)) in enumerate(data.arc))

    c13 = ExaModels.constraint!(w, c9, g.bus => -pg[i] for (i, g) in enumerate(data.gen))
    c14 = ExaModels.constraint!(w, c10, g.bus => -qg[i] for (i, g) in enumerate(data.gen))

    return ExaModels.ExaModel(w; prod = true),
    (va, vm, pg, qg, p, q),
    (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)

end

function _jump_ac_power_struct_model(backend, filename)

    data = parse_ac_power_struct_data(filename, backend)

    model = JuMP.Model()
    #JuMP.set_optimizer_attribute(model, "print_level", 0)

    JuMP.@variable(model, va[i in 1:length(data.bus)])
    JuMP.@variable(
        model,
        data.vmin[i] <= vm[i in 1:length(data.bus)] <= data.vmax[i],
        start = 1.0
    )

    JuMP.@variable(
        model,
        data.pmin[i] <= pg[i in 1:length(data.gen)] <= data.pmax[i],
    )
    JuMP.@variable(
        model,
        data.qmin[i] <= qg[i in 1:length(data.gen)] <= data.qmax[i],
    )

    JuMP.@variable(
        model,
        -data.rate_a[i] <= p[i in 1:length(data.arc)] <= data.rate_a[i]
    )
    JuMP.@variable(
        model,
        -data.rate_a[i] <= q[i in 1:length(data.arc)] <= data.rate_a[i]
    )

    JuMP.@NLobjective(
        model,
        Min,
        sum(
            g.c[1] * pg[i]^2 + g.c[2] * pg[i] + g.c[3]
            for (i, g) in enumerate(data.gen)
        )
    )

    c1 = []
    c2 = []
    c3 = []
    c4 = []
    c5 = []
    c6 = []
    c7 = []
    c8 = []
    c9 = []
    c10 = []

    for b in data.ref_buses
        push!(c1, JuMP.@NLconstraint(model, va[b] == 0))
    end

    # Branch power flow physics and limit constraints
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(
            c2,
            JuMP.@NLconstraint(
                model,
                p[i.fidx] == b.c5 * vm[b.fbus]^2 +
                b.c3 * (vm[b.fbus] * vm[b.tbus] * cos(va[b.fbus] - va[b.tbus])) +
                b.c4 * (vm[b.fbus] * vm[b.tbus] * sin(va[b.fbus] - va[b.tbus]))
            )
        )
    end
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(
            c3,
            JuMP.@NLconstraint(
                model,
                q[i.fidx] ==
                -b.c6 * vm[b.fbus]^2 -
                b.c4 * (vm[b.fbus] * vm[b.tbus] * cos(va[b.fbus] - va[b.tbus])) +
                b.c3 * (vm[b.fbus] * vm[b.tbus] * sin(va[b.fbus] - va[b.tbus]))
            )
        )
    end
    # To side of the branch flow
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(
            c4,
            JuMP.@NLconstraint(
                model,
                p[i.tidx] ==
                b.c7 * vm[b.tbus]^2 +
                b.c1 * (vm[b.tbus] * vm[b.fbus] * cos(va[b.tbus] - va[b.fbus])) +
                b.c2 * (vm[b.tbus] * vm[b.fbus] * sin(va[b.tbus] - va[b.fbus]))
            )
        )
    end
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(
            c5,
            JuMP.@NLconstraint(
                model,
                q[i.tidx] ==
                -b.c8 * vm[b.tbus]^2 -
                b.c2 * (vm[b.tbus] * vm[b.fbus] * cos(va[b.tbus] - va[b.fbus])) +
                b.c1 * (vm[b.tbus] * vm[b.fbus] * sin(va[b.tbus] - va[b.fbus]))
            )
        )
    end
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(
            c6,
            JuMP.@NLconstraint(
                model,
                b.angmin <= va[b.fbus] - va[b.tbus] <= b.angmax
            )
        )
    end
    # Apparent power limit, from side and to side
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(c7, JuMP.@NLconstraint(model, p[i.fidx]^2 + q[i.fidx]^2 <= b.ratea^2))
    end
    for (i, b) in zip(data.arc_idxs, data.branch)
        push!(c8, JuMP.@NLconstraint(model, p[i.tidx]^2 + q[i.tidx]^2 <= b.ratea^2))
    end

    bus_arcs = [[] for i in 1:length(data.bus)]
    for (i, (l, a_bus, j)) in enumerate(data.arc)
        push!(bus_arcs[a_bus], i)
    end
    bus_gens = [[] for i in 1:length(data.bus)]
    for (i, g) in enumerate(data.gen)
        push!(bus_gens[g.bus], i)
    end
    for (i, b) in enumerate(data.bus)
        push!(
            c9,
            JuMP.@NLconstraint(
                model,
                sum(p[a] for a in bus_arcs[i]) ==
                sum(pg[g] for g in bus_gens[i]) -
                b.pd -
                b.gs * vm[i]^2
            )
        )
    end
    for (i, b) in enumerate(data.bus)
        push!(
            c10,
            JuMP.@NLconstraint(
                model,
                sum(q[a] for a in bus_arcs[i]) ==
                sum(qg[g] for g in bus_gens[i]) -
                b.qd +
                b.bs * vm[i]^2
            )
        )
    end

    return model,
    (va, vm, pg, qg, p, q),
    (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
end
