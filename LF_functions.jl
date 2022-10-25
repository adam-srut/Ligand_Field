function S_lm(l::Int, m::Int, θ, ϕ)
	#= Computes a structure factor Yₗᵐ(θ,ϕ) for a ligand at position [θ,ϕ]=#
    S = 4*pi/(2*l+1)*SphericalHarmonics.sphericalharmonic(θ, ϕ, l=l, m=m)
    return S
end

function gaunt(l3::Int, m3::Int, m1::Int, m2::Int)
    #= Computes the Gaunt integral in a form:
        < Y₂,ₘ₁ | Yₗ₃,ₘ₃ | Y₂,ₘ₂ > =#
    G = (-1)^(m1)*sqrt((25*(2*l3+1))/(4*pi)) * 
        wigner3j(Float64, 2, 2, l3, 0, 0, 0) * 
        wigner3j(Float64, 2, 2, l3, -m1, m2, m3)
    return G
end

function construct_VLF(m1::Int, m2::Int, ligands_θ::Array, ligands_ϕ::Array)
    #= Function returns factors for Slater-Condon parameters
        for given element of the total Hamiltonian < 2,m1 | V^LF | 2,m2 >.
    =#
    NL = length(ligands_θ)
    V_LF = 0
    (F0, F2, F4) = (0, 0, 0)
    for i_lig in 1:NL
        F0 += S_lm(0, 0, ligands_θ[i_lig], ligands_ϕ[i_lig]) * gaunt(0, 0, m1, m2)
        for m in -2:2
            F2 += S_lm(2, m, ligands_θ[i_lig], ligands_ϕ[i_lig]) * gaunt(2, m, m1, m2)
        end
        for m in -4:4
            F4 += S_lm(4, m, ligands_θ[i_lig], ligands_ϕ[i_lig]) * gaunt(4, m, m1, m2)
        end
    end
    return (F0,F2,F4)
end
            
function spherical_to_cartesian(θ::Number, ϕ::Number)
    x = cos(ϕ)*sin(θ)
    y = sin(ϕ)*sin(θ)
    z = cos(θ)
    return(x,y,z)
end

function plot_ligs(ligands_θ::Array, ligands_ϕ::Array)
    xyz_ligands = [ [spherical_to_cartesian(θ, ϕ)...] for (θ, ϕ) in zip(ligands_θ, ligands_ϕ) ]
    xyz_ligands = reduce(hcat, xyz_ligands)'
    xyz_lines = Matrix{Float64}(undef, length(eachrow(xyz_ligands))*2, 3)
    for (i,coord) in enumerate(eachrow(xyz_ligands))
        xyz_lines[2*i-1,:] = coord
        xyz_lines[2*i,:] = [0,0,0]
    end
    plot(xyz_lines[:,1], xyz_lines[:,2], xyz_lines[:,3],
        color = :black, legend = false)
    scatter!(xyz_ligands[:,1], xyz_ligands[:,2], xyz_ligands[:,3],m = (12,0.99, :blues, Plots.stroke(0.5)))
    scatter!([0],[0],[0], m = (10,20,:orange, Plots.stroke(3)))
    plot!(showaxis=false, grid=false, ticks=false, zlims=(-1,1), ylims=(-1,1), xlims=(-1,1))
end

function construct_Htot_DF(ligands_θ::Array, ligands_ϕ::Array)
    V_LF = Matrix{String}(undef,5,5)
    for i in -2:2
        for j in -2:2
            V_LF[i+3, j+3] = ""
            factors = real.(construct_VLF(i, j, ligands_θ, ligands_ϕ))
            for f in 1:3
                if abs(factors[f]) < 1e-10
                    continue
                else
                    factors[f] > 0 ? V_LF[i+3, j+3] *= "+" : V_LF[i+3, j+3] *= "-"
                    V_LF[i+3, j+3] *= @sprintf "%.2f F%d" abs(factors[f]) (f-1)*2
                end
            end
        end
    end
    table = DataFrame( 
        "d₋₂" => V_LF[:,1], 
        "d₋₁" => V_LF[:,2],
        "d₀" => V_LF[:,3],
        "d₁" => V_LF[:,4],
        "d₂" => V_LF[:,5])
    return table
end

function compute_Htot(ligands_θ::Array, ligands_ϕ::Array;
        f0::Number, f2::Number, f4::Number)
    V_LF = zeros(5,5)
    for i in -2:2
        for j in -2:2
            factors = real.(construct_VLF(i, j, ligands_θ, ligands_ϕ))
            for (i_f,f) in enumerate([f0,f2,f4])
                if abs(factors[i_f]) < 1e-20
                    continue
                else
                    V_LF[i+3, j+3] += factors[i_f]*f
                end
            end
        end
    end
    return V_LF
end

function draw_splitting(eig_vals::Array)
    Drawing(600,400, "d-orb_splitting.png")
    origin()
    setcolor("gray")
    setline(0.8)
    line(Point(-400,0),Point(400,0), :stroke)
    sethue("black")
    setline(2.5)
    fontsize(14)
    for (shift, orb) in enumerate(["d₋₂", "d₋₁", "d₀", "d₁", "d₂"])
        Luxor.text(orb, Point(30*shift-300, 25),halign=:center, valign=:center)
    end
    for i in 1:5
        line(Point(30*i-310,0),Point(30*i-290,0), :stroke)
    end
    for (i,val) in enumerate(eig_vals)
        shift = 0
        if i > 1
            if (val - eig_vals[i-1]) < 0.1
                shift = 50
                if i > 2
                    (val - eig_vals[i-2]) < 0.1 ? shift = 100 : shift = 50
                end
            else 
                shift = 0
            end
        end
        line(Point(-50+shift,-val*40), Point(-10+shift,-val*40), :stroke)
    end
    Luxor.arrow(Point(-120,150),Point(-120,-150), arrowheadlength=22, linewidth=2.2)
    fontsize(20)
    Luxor.text("E", Point(-150,-120))
    finish()
    preview()
end

function do_splitting(ligands_θ::Array, ligands_ϕ::Array;
        f0::Number, f2::Number, f4::Number)
    Htot = compute_Htot(ligands_θ, ligands_ϕ, f0=f0, f2=f2, f4=f4)
    (vals, vecs) = eigen(Htot)
    draw_splitting(vals)
end

function construct_Htot(ligands_θ::Array, ligands_ϕ::Array)
    V_LF = "\\begin{smallmatrix}"#  & d_{-2} & d_{-1} & d_0 & d_1 & d_2 \\\\"
    for i in -2:2
        #V_LF *= "d_{$i} &"
        for j in -2:2
            factors = real.(construct_VLF(i, j, ligands_θ, ligands_ϕ))
            for f in 1:3
                if abs(factors[f]) < 1e-10
                    continue
                else
                    factors[f] > 0 ? V_LF *= " + " : V_LF *= " - "
                    V_LF *= @sprintf "%.2f" abs(factors[f])
                    num = Int((f-1)*2)
                    V_LF *= "\\cdot F_$(num)"
                end
            end
            V_LF *= "&"
        end
        V_LF *= "\\\\"
    end
    V_LF *= "\\end{smallmatrix}"
    return latexstring(V_LF)
end

function print_ExE(directory::String)
    cmd = "grep -A 15 'ABSORPTION SPECTRUM' single-point-calc.out | tail -n 15"
    run(Cmd(`/bin/bash -c $cmd`, dir=directory))
end

function print_OrbE(directory::String)
    cmd = "grep -A 3 -B 1 'ORBITAL ENERGIES' single-point-calc.out"
    run(Cmd(`/bin/bash -c $cmd`, dir=directory))
    cmd = "grep -A 47 'ORBITAL ENERGIES' single-point-calc.out | tail -n 5"
    run(Cmd(`/bin/bash -c $cmd`, dir=directory))
end

function run_calculation(dirname::String)
	cmd = "\$ORCADIR/orca single-point-calc.inp > single-point-calc.out"
    run(Cmd(`/bin/bash -c $cmd`, dir=dirname));
end
