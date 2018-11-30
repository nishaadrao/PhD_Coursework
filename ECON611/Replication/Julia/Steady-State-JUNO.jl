using QuantEcon, Interpolations, Plots, LaTeXStrings, DataFrames
import PyPlot

######################################################################
# CALIBRATION
######################################################################

const β     = 0.986;
const Rbar  = 1.005;
const γ     = 2;
const μ     = 1.2;
const ψ     = 2;
const Wbar  = 1/μ; #steady state wage

# Tax function (only rich pay tax)
function tax(x)
    if x==z_grid[3]
        return 1
    else return 0
    end
end

# Markov chain for z
MC = rouwenhorst(3, 0.96566, sqrt(0.01695/(1-0.96566^2)), 0);

# z grid
z_grid  = [0.49229735012925,1,2.03129267248230];

# bond grid
b_min  = 0;
b_max  = 75;
b_grid = collect(linspace(b_min,b_max,200));


######################################################################
# EGM function
######################################################################

function egm{TR<:Real,TF<:AbstractFloat}(G::Array{TF,2},
             b_grid::Array{TF,1},
             z_grid::Array{TF,1},
             P::Array{Float64,2},
             β::TR,
             γ::TR,
             ψ::TR,
             R::TR,
             W::TR,
             D::TR,
             τ::TR)

    # Allocate memory for value of today's consumption on endogenous grid points,...
    # and the endogenous grid itself,...
    # and the updated guess of the policy function for each state pair.
    c      = zeros(length(z_grid),length(b_grid))
    b_eg   = zeros(length(z_grid),length(b_grid))
    KG     = zeros(length(z_grid),length(b_grid))

    # Compute today's optimal consumption using Euler equation
    for (i, z) in enumerate(z_grid)
        for (j,b) in enumerate(b_grid)
            c[i,j] = (β*R*P[i,:]'*(G[:,j].^(-γ)))^(-1/γ)
        end
    end

    # Compute endogenous bond grid using the budget constraint
    for (i,z) in enumerate(z_grid)
        for(j,b) in enumerate(b_grid)
            b_eg[i,j] = c[i,j] + b_grid[j]/R - c[i,j]^(-γ/ψ)*(W*z_grid[i])^(1+1/ψ) - D + τ*tax(z_grid[i])
        end
    end

    # Compute today's optimal consumption for agents with binding constraint
    for (i,z) in enumerate(z_grid)
        for(j,b) in enumerate(b_grid)

            if b_grid[j] <= b_eg[i,1]
                c[i,j] = (b_grid[j]+D-τ*tax(z_grid[i])+sqrt((b_grid[j]+D-τ*tax(z_grid[i]))^2+4*(W*z_grid[i])^(1+1/ψ)))/2
            else
                c[i,j] = c[i,j]
            end
        end
    end


    # Interpolate the updated guess [NEED TO CHANGE THIS TO CUBIC SPLINE!]
    for (i,z) in enumerate(z_grid)
        KG_f      = LinInterp(b_eg[i,:],c[i,:])
        KG[i,:]   = KG_f.(b_grid)'
    end

    return KG

end

######################################################################
# Initial guess of policy function
######################################################################

function G_init{TR<:Real}(W::Real,D::TR,τ::TR)

    G0   = zeros(length(z_grid),length(b_grid))

    for (i,z) in enumerate(z_grid)
            for(j,b) in enumerate(b_grid)
                G0[i,j] = (b_grid[j]+D-τ*tax(z_grid[i])+sqrt((b_grid[j]+D-τ*tax(z_grid[i]))^2+4*(W*z_grid[i])^(1+1/ψ)))/2
            end
    end
    return G0

end

######################################################################
# Compute policy function
######################################################################

function get_policy{TR<:Real}(G0,R::TR,W::TR, D::TR, τ::TR, tol=1e-8, maxiter=1000, err=1, i=0)
    G       = G0    # Initial guess

    while i<=maxiter && err > tol

        Gnew = egm(G,b_grid,z_grid,MC.p,β,γ,ψ,R,W,D,τ)

        err  = maximum(abs, Gnew - G)
        G    = Gnew

        i    = i+1
    end
    return G

end


######################################################################
# Simulate bond distribution
######################################################################

function mc_sample_path(P, init=2, sample_size=1000)
    X = Array{Int64}(sample_size) # allocate memory
    Z = zeros(length(X))

    X[1] = init

    # === convert each row of P into a distribution === #
    n = size(P)[1]
    P_dist = [DiscreteRV(vec(P[i,:])) for i in 1:n]

    # === generate the sample path for the state (i.e. X \in {1,2,3}) === #
    for t in 1:(sample_size - 1)
        X[t+1] = rand(P_dist[X[t]])
    end

    # === get the z associated with each value of the state === #
    for t in 1:(sample_size)
        if X[t]==1
            Z[t]=z_grid[1]
        elseif X[t]==2
            Z[t]=z_grid[2]
        else
            Z[t]=z_grid[3]
        end
    end

    return Z
end

# simulate a long vector of states
Z=mc_sample_path(MC.p,2,1000000)

# Function to get associated bond holdings, using policy function
function get_bonds{TR<:Real}(Z,G::AbstractArray,R::Real,W::Real,D::TR,τ::TR,c1,c2,c3,init=b_grid[1])
    B = zeros(length(Z)) #allocate memory

    # This is what c1,c2,c3 will be!
    #c1 = LinInterp(b_grid,G[1,:]);
    #c2 = LinInterp(b_grid,G[2,:]);
    #c3 = LinInterp(b_grid,G[3,:]);

    B[1] = init

    for t in 1:(length(B)-1)

        if Z[t]==z_grid[1]
            B[t+1]=R*(B[t]+(W*Z[t])^(1+1/ψ)*c1(B[t])^(-γ/ψ)-c1(B[t])+D-τ*tax(z_grid[1]))
        elseif Z[t]==z_grid[2]
            B[t+1]=R*(B[t]+(W*Z[t])^(1+1/ψ)*c2(B[t])^(-γ/ψ)-c2(B[t])+D-τ*tax(z_grid[2]))
        else
            B[t+1]=R*(B[t]+(W*Z[t])^(1+1/ψ)*c3(B[t])^(-γ/ψ)-c3(B[t])+D-τ*tax(z_grid[3]))
        end
    end

    return B
end

# Get distribution of consumption
function get_C(B,c1,c2,c3)

    C = zeros(length(B)) #allocate memory

    for t in 1:(length(C))

            if Z[t]==z_grid[1]
                C[t]=c1(B[t])
            elseif Z[t]==z_grid[2]
                C[t]=c2(B[t])
            else
                C[t]=c3(B[t])
            end
    end

    return C
end

######################################################################
# Solve for steady state
######################################################################
#function solve_ss(D=0,tol=1e-8,maxiter=1000,err=1,i=0)

W       = Wbar   # Initial guess
R       = Rbar
D       = 0
τ       = 0
tol     = 1e-8
maxiter = 1000
err     = 1
i       = 0


    while i<=maxiter && err > tol

        global G, W, R, C, D, B, τ       # Ensures the loop spits out these things (cf. local)

        X    = [R,W,D,τ]

        G0   = G_init(W,D,τ)             # Compute initial guess of policy fn

        G    = get_policy(G0,R,W,D,τ)    # Compute policy function

        c1 = LinInterp(b_grid,G[1,:])    # optimal policy *function* for low productivity
        c2 = LinInterp(b_grid,G[2,:])
        c3 = LinInterp(b_grid,G[3,:])

        B  = get_bonds(Z,G,R,W,D,τ,c1,c2,c3) # Simulate distribution of bonds

        C  = get_C(B,c1,c2,c3)           # Get consumption distribution

        # CHECK LABOR MARKET CLEARING AND ADJUST WAGES
        L    = C.^(-γ/ψ).*((W*Z).^(1/ψ)); # Compute lsupply using intratemporal sub condition
        zL   = L.*Z;
        Lbar = mean(zL)

        Wtil = W + 0.8*(mean(C)-Lbar)

        Wnew = 0.5*Wtil + 0.5*W

        # CHECK BOND MARKET CLEARING AND ADJUST INTEREST RATE
        Rtil = R - 0.001*(mean(B)-1.4*4*mean(C))

        Rnew = 0.5*Rtil + 0.5*R

        Dnew = mean(C)*(1-W)             # Compute prelim updated guess of divs

        τnew = 1.4*4*mean(C)*(Rnew-1)/Rnew     # Compute prelim updated guess of τ

        Xnew = [Rnew,Wnew,Dnew,τnew]

        err  = maximum(abs, Xnew - X)  # Compute error

        R    = Rnew
        W    = Wnew
        D    = 0.5*Dnew + 0.5*D        # Update guess of dividends
        τ    = 0.5*τnew + 0.5*τ         # Update guess of τ

        i    = i+1
    end
#end
