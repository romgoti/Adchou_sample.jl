

"""
solve_2states(T,Aprod,2,γ,α,δ,ρ,K,z1,λ1,λ2,Δ ,print_ind,amin,amax,maxit,crit)

Solves for the 2-states equilibrium for given parameter values and options for estimation
(V, K, r, g, cf, sf , da, w, L) 
"""

function solve_2states(
    T=250,        # number of a points 
    Aprod = .1,  # Aggregate productivity
    z2 = 2,      # multiple of low state income z1 => z2 = z1 * z2
    γ = 1.2,       # CRRA utility with parameter γ
    α = 1/3,   # Production function F = K^alpha * L^(1-alpha) 
    δ = 0.0,  # Capital depreciation => for delta=0 put guess Kss ≈ 13
    ρ = 0.05,   # discount rate
    K = 12.0,      # initial aggregate capital. It is important to guess a value close to the solution for the algorithm to converge (otherwise returns an error)
    z1 = 1,     # low-state z
    λ1 = 1/4,   # transition rates
    λ2 = 1/4,
    Δ = 1000,    # time increment 
    print_ind = 0,     # option to print convergence status
    amin = 0,    # borrowing constraint
    amax = 8.,    # highest a in grid
    maxit = 100,    # max nb of iterations (inner loop)
    crit = 10^(-9.)  # tolerance criterion HJB loop
)

    z = [z1,z2]
    λ = [λ1,λ2]
    z2 = z2*z1;
    zbar = (z1*λ1 + z2*λ2)/(λ1 + λ2); # average individual productivity
    a = range(amin,amax,T);  #wealth vector
    da = (amax-amin)/(T-1);    # wealth increment
    aa = a*ones(1,2);
    zz = ones(T,1)*z';

    # utility function
    u(x) = x^(1-γ)/(1-γ)

    #Finite difference approximation of the partial derivatives
    # initialize the vectors 
    dVf = zeros(T,2)
    dVb = zeros(T,2)
    v0 = zeros(T,2)
    c = zeros(T,2)
    cf = zeros(T,2)
    Id = sparse(1:T,1:T,ones(T))
    ID = sparse(1:T*2,1:T*2,ones(T*2))

    # Build the transition matrix in continuous time
    Lswitch = spdiagm(0 => -[fill(λ1,(T,1))...,fill(λ2,(T,1))...],-T => ones(T)λ2,T => ones(T)λ1)
    L = zeros(2T,2T)
    LT = zeros(2T,2T)
    dist = zeros(maxit) 
    g = zeros(T,2)

    
    Ir = 100;       # max number of iterations to converge (outer loop)
    crit_S = 10^(-8) # inner loop criterion

    # initial guesses and initialize objects
    r = 0.05;
    w = 0.05;
    rmin = 0.01;
    rmax = 0.99*ρ;
    r_list = zeros(Ir)
    rmin_list = zeros(Ir)
    rmax_list = zeros(Ir)
    KD = zeros(Ir)
    KS = zeros(Ir)
    V = zeros(T,2)
    sf = zeros(T,2)
    V_r = zeros(T,2)
    id = 0
    i_fix = 0

    #----------------------------------------------------
    #INITIAL GUESS
    v0[:,1] = u.(z[1] .+ r*a)/ρ;
    v0[:,2] = u.(z[2] .+ r*a)/ρ;

    #-----------------------------------------------------
    #OUTER LOOP
    for ir=1:Ir

        if print_ind == 1
            println("Main loop for ",ir);
        end
        id = ir

        # store past r guesses
        r_list[ir] = r;
        rmin_list[ir] = rmin;
        rmax_list[ir] = rmax;

        # compute K demand and associated wage
        KD[ir] = ((α*Aprod)/(r+δ))^(1/(1-α))*zbar
        w = (1-α)*Aprod * KD[ir]^(α)*zbar^(-α);  

        if ir > 1
            v0 = V_r
        end

        v=v0

        # HAMILTON-JACOBI-BELLMAN EQUATION %
        # INNER LOOP

        # inefficient => just to solve for the partial equilibrium once if needed
        # v,L = solve_forw_dif(v,da, T, w, r, amax, amin, z, zz, aa, γ, u, Lswitch,Δ,ρ,print_ind) 

        for n=1:maxit
            V = v;
            #V_n
            
            # forward difference
            dVf[1:T-1,:] = (V[2:T,:]-V[1:T-1,:])/da;
            dVf[T,:] = (w*z .+ r*amax).^(-γ); #will never be used, but impose state constraint a<=amax just in case
            # backward difference
            dVb[2:T,:] = (V[2:T,:]-V[1:T-1,:])/da;
            dVb[1,:] = (w*z .+ r*amin).^(-γ);  #state constraint boundary condition
    
            #I_concave = Vab > Vaf; indicator whether value function is concave (problems arise if this is not the case)
    
            #consumption and savings with forward difference
            cf = dVf.^(-1/γ);
            sf = w*zz + r*aa - cf;
    
            #consumption and savings with backward difference
            cb = dVb.^(-1/γ);
            sb = w*zz + r*aa - cb;
            #consumption and derivative of value function at steady state
            c0 = w*zz + r*aa;
            Va0 = c0.^(-γ);
    
            # dV_upwind makes a choice of forward or backward differences based on
            # the sign of the drift
            If = sf .> 0; #positive drift --> forward difference
            Ib = sb .< 0; #negative drift --> backward difference
    
            I0 = (ones(T,2)-If-Ib); #at steady state
    
            #make sure backward difference is used at amax
            #     Ib[I,:] = 1; If[I,:] = 0;
            #STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS sf > 0:
            #already taken care of automatically
    
            Va_Upwind = dVf.*If + dVb.*Ib + Va0.*I0; #important to include third term
    
            cf = Va_Upwind.^(-1/γ);
            uc = u.(cf);
    

            # SOLVE ρv = u(v) + Av (in paper)

            #CONSTRUCT MATRIX A
            X = - min.(sb,0)/da;
            Y = - max.(sf,0)/da + min.(sb,0)/da;
            Z = max.(sf,0)/da;
    
            # transition matrix including policy response
            L = spdiagm(0 => [Y[:,1]...,Y[:,2]...], -1 => [X[2:T,1]...,0,X[2:T,2]...], 1 => [Z[1:T-1,1]...,0,Z[1:T-1,2]...]);
            L += Lswitch;
            
            if maximum(abs.(sum(L,dims=2))) > 10^(-9)
                print("Improper Transition Matrix")
                break
            end
    
            B = (1/Δ + ρ)*ID - L;
    
            u_stacked = reshape(uc,2T,1);
            V_stacked = reshape(V,2T,1);
    
            BLAS.axpy!(2T,1/Δ,V_stacked,1,u_stacked,1);
    
            V_stacked = B\u_stacked[:,1];#B\b[:,1]; #SOLVE SYSTEM OF EQUATIONS
    
            V = reshape(V_stacked,T,2);
    
            Vchange = V - v;
            v = V;
    
            dist[n] = maximum(abs.(Vchange));
            if dist[n]<crit
                if print_ind == 1
                    println("Value Function Converged, Iteration = ",n);
                end
                break
            end
        end

        # FOKKER-PLANCK EQUATION %
        LT = L';
        b = zeros(2T,1)

        #need to fix one value, otherwise matrix is singular
        i_fix = floor(Int64,T/2);
        b[i_fix]=.1;

        for j=1:T*2
            LT[i_fix,j]=0;
        end

        LT[i_fix,i_fix]=1;

        #Solve linear system : 0 = A'g
        gg = LT\b;
        g_sum = gg'*ones(2T,1)*da;
        gg = gg./g_sum;

        g = reshape(gg,T,2);

        # check1 = g[:,1]'*ones(T,1)*da;
        # check2 = g[:,2]'*ones(T,1)*da;
        
        # Update aggregate capital
        g_r = g;
        KS[ir] = sum(g.*a*da); # K supply
        V_r = V;
        S = KS[ir] - KD[ir] # excess demand or excess supply of capital
        
        if print_ind == 1
            println(S)
        end

    
        if S>crit_S 
            if print_ind == 1
                print("Excess Supply")
            end
            rmax = r 
            r = .5(r + rmin)
        elseif S<-crit_S
            if print_ind == 1
                print("Excess Demand")
            end
            rmin = r 
            r = .5(r + rmax)
        elseif abs(S) <= crit_S
            if print_ind == 1
                print("Equilibrium at ")
                print(r) 
            end

            break
        end
        
    end

    K = KS[id] # final K scalar

    # smooth L where needed
    # L[i_fix,i_fix] = (L[i_fix-1,i_fix-1]+L[i_fix+1,i_fix+1])/2
    # L[i_fix+1,i_fix] = (L[i_fix,i_fix-1]+L[i_fix+2,i_fix+1])/2
    # L[T+i_fix,i_fix] = (L[T+i_fix-1,i_fix-1]+L[T+i_fix+1,i_fix+1])/2

    return OrderedDict(:V => V , :K  => K , :r => r, :g => g, :cf => cf, :sf => sf ,:da => da, :w => w, :L => L) 
end


# function to solve the partial equilibrium (inner loop) separately, can be useful if needed to compare interest rates 
function solve_forw_dif(v,da, T, w, r, amax, amin, z, zz, aa, γ, u, Lswitch,Δ,ρ,print_ind,maxit,crit)
    
    dVf = zeros(T,2);
    dVb = zeros(T,2);
    cf = zeros(T,2);
    sf = zeros(T,2);
    ID = sparse(1:T*2,1:T*2,ones(T*2));
    dist = zeros(maxit);
    L = zeros(2T,2T);

    for n=1:maxit
        V = v;
        
        # forward difference
        dVf[1:T-1,:] = (V[2:T,:]-V[1:T-1,:])/da;
        dVf[T,:] = (w*z .+ r*amax).^(-γ); #will never be used, but impose state constraint a<=amax just in case
        # backward difference
        dVb[2:T,:] = (V[2:T,:]-V[1:T-1,:])/da;
        dVb[1,:] = (w*z .+ r*amin).^(-γ);  #state constraint boundary condition

        #consumption and savings with forward difference
        cf = dVf.^(-1/γ);
        sf = w*zz + r*aa - cf;

        #consumption and savings with backward difference
        cb = dVb.^(-1/γ);
        sb = w*zz + r*aa - cb;
        #consumption and derivative of value function at steady state
        c0 = w*zz + r*aa;
        Va0 = c0.^(-γ);

        If = sf .> 0; #positive drift --> forward difference
        Ib = sb .< 0; #negative drift --> backward difference

        I0 = (ones(T,2)-If-Ib); #at steady state

        Va_Upwind = dVf.*If + dVb.*Ib + Va0.*I0; #important to include third term

        cf = Va_Upwind.^(-1/γ);
        uc = u.(cf);

        #CONSTRUCT MATRIX A
        X = - min.(sb,0)/da;
        Y = - max.(sf,0)/da + min.(sb,0)/da;
        Z = max.(sf,0)/da;

        L = spdiagm(0 => [Y[:,1]...,Y[:,2]...], -1 => [X[2:T,1]...,0,X[2:T,2]...], 1 => [Z[1:T-1,1]...,0,Z[1:T-1,2]...]);
        L += Lswitch;
        
        if maximum(abs.(sum(L,dims=2))) > 10^(-9)
            print("Improper Transition Matrix")
            break
        end

        B = (1/Δ + ρ)*ID - L;

        u_stacked = reshape(uc,2T,1);
        V_stacked = reshape(V,2T,1);

        BLAS.axpy!(2T,1/Δ,V_stacked,1,u_stacked,1);

        V_stacked = B\u_stacked[:,1];#B\b[:,1]; #SOLVE SYSTEM OF EQUATIONS

        V = reshape(V_stacked,T,2);

        Vchange = V - v;
        v = V;

        dist[n] = maximum(abs.(Vchange));
        if dist[n]<crit
            if print_ind == 1
                println("Value Function Converged, Iteration = ",n);
            end
            break
        end
    end

    return OrderedDict(:v => v, :L => L, :c => cf, :s => sf)
end


# HANDLING NON-CONVEXITIES SECTION 6
# solves the housing model
function solve_housing(
    T=500,         # number of a points 
    scenario = 1,  # starting point for distribution of assets
    # since this model has a "poverty trap" the resulting distribution of asset holdings depends on how it is simulated 
    # scenario 1 : SS with some workers constrained far under the down-payment for housing good and some well-above
    # scenario 2 : SS with all workers constrained under the down-payment for housing 
    # scenario 3 : SS with all workers holding some housing

    print_ind = 0,   # print status estimation
    s = 2,        # CRRA parameter
    r = 0.035,    # scalar for bond interest rate
    α = 1/3,
    η = 0.2,    # housing weight in the utility
    p = 1,      # price of housing
    ϕ = 2,          # ϕ = 1/(1-θ) with θ the down payment share (here 1/2)
    hmin = 2.3,     # minimum size of housing good
    ρ = 0.05,     
    z1 = .1,
    z2 = .135,
    λ1 = 0.5,
    λ2 = 0.5,
    amin = 0,
    amax = 3
)
        
    astar = p*hmin/ϕ;          # non-linear point above which HHs can make the down payment
    λ = [λ1,λ2];
    z = [z1,z2];

    a = range(amin,amax,T)';
    da = (amax-amin)/(T-1);
    aa = [a ; a]';
    zz = ones(T,1)*z';

    maxit= 120;
    crit = 10^(-6);
    Delta = 500;

    # initialize
    dVf = zeros(T,2);
    dVb = zeros(T,2);
    c = zeros(T,2);
    A = zeros(2T,2T) 

    # transition matrix in continuous time
    Aswitch = spdiagm(0 => [-λ1*ones(T) ; -λ2*ones(T)],-T => λ2*ones(T),T => λ1*ones(T));

    # housing choice function => wealth allocation
    h = min.((α*η/(r*p))^(1/(1-α)) + hmin,ϕ*aa/p);
    h = h.*(h.>=hmin);
    f = η*(max.(h.-hmin,0)).^α - r*p*h;

    # initial guess
    v0 = (zz + r*aa).^(1-s)/(1-s)/ρ;
    v = v0;

    for n = 1:maxit

        V = v;
        # V_n(:,:,n)=V;

        # % forward difference
        dVf[1:T-1,:] = (V[2:T,:]-V[1:T-1,:])/da;
        dVf[T,:] = (z + f[T,:] .+ r*amax).^(-s); #%will never be used, but impose state constraint a<=amax just in case
        
        # % backward difference
        dVb[2:T,:] = (V[2:T,:]-V[1:T-1,:])/da;
        dVb[1,:] = (z + f[1,:] .+ r*amin).^(-s); #%state constraint boundary condition
        
        # %consumption and savings with forward difference
        cf = (max.(dVf,10^(-10))).^(-1/s);
        ssf = zz + f + r.*aa - cf;
        Hf = cf.^(1-s)/(1-s) + dVf.*ssf;
        
        # %consumption and savings with backward difference
        cb = (max.(dVb,10^(-10))).^(-1/s);
        ssb = zz + f + r.*aa - cb;
        Hb = cb.^(1-s)/(1-s) + dVb.*ssb;
        
        # %consumption and derivative of value function at steady state
        c0 = zz + f + r.*aa; 
        
    # %     % dV_upwind makes a choice of forward or backward differences based on    
        Ineither = (1 .-(ssf.>0)) .* (1 .-(ssb.<0));
        Iunique = (ssb.<0).*(1 .-(ssf.>0)) + (1 .-(ssb.<0)).*(ssf.>0);
        Iboth = (ssb.<0).*(ssf.>0);
        Ib = Iunique.*(ssb.<0) + Iboth.*(Hb.>=Hf);
        If = Iunique.*(ssf.>0) + Iboth.*(Hf.>=Hb);
        I0 = Ineither;
        
        c = cf.*If + cb.*Ib + c0.*I0;
        uu = c.^(1-s)/(1-s);
        
        # %CONSTRUCT MATRIX
        X = - Ib.*ssb/da;
        Y = - If.*ssf/da + Ib.*ssb/da;
        Z = If.*ssf/da;
        
        # A1=spdiags(Y(:,1),0,T,T)+spdiags(X(2:T,1),-1,T,T)+spdiags([0;Z(1:T-1,1)],1,T,T);
        # A2=spdiags(Y(:,2),0,T,T)+spdiags(X(2:T,2),-1,T,T)+spdiags([0;Z(1:T-1,2)],1,T,T);
        A = spdiagm(0 => [Y[:,1] ; Y[:,2]],-1 => [X[2:T,1] ; 0 ; X[2:T,2]],1 => [Z[1:T-1,1] ; 0 ; Z[1:T-1,2]]);
        A = A + Aswitch;
        
        if maximum(abs.(sum(A,dims = 2)))>10^(-12)
            println("Improper Transition Matrix")
            break
        end
        
        B = (1/Delta + ρ)*Diagonal(ones(2*T)) - A;
        
        u_stacked = [uu[:,1];uu[:,2]];
        V_stacked = [V[:,1];V[:,2]];
        
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b; #%SOLVE SYSTEM OF EQUATIONS
        
        V = [V_stacked[1:T] V_stacked[T+1:2*T]];
        
        Vchange = V - v;
        v = V;

        dist = maximum(abs.(Vchange));
        if dist < crit
            if print_ind == 1
                println("Value Function Converged, Iteration = ",n)               
            end
            break
        end
    end


    # %%%%%%%%%%%%%%%%%%%%%%%%%%
    # % FOKKER-PLANCK EQUATION %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%
    AT = A';

    if scenario == 1
         # %INITIAL DISTRIBUTION 1:
         gg0 = ones(2*T,1);       
    elseif scenario == 2
        # %INITIAL DISTRIBUTION 2
        j = Int(T/2); gg0 = [zeros(j,1);ones(T-j,1);zeros(j,1);ones(T-j,1)];
    else scenario == 3
        # %INITIAL DISTRIBUTION 3
        j = 10; gg0 = [ones(j,1);zeros(T-j,1);ones(j,1);zeros(T-j,1)];
    end

    g_sum = gg0'*ones(2*T,1)*da;
    gg0 = gg0./g_sum;

    gg=gg0;
    N=1000; dt=10;
    for n=1:N
        # %Implicit method in Updating Distribution.
        gg0 = gg;
        gg= (Diagonal(ones(2*T)) - AT*dt)\gg;
        g_dist=maximum(abs.(gg-gg0));
    end

    adot = zz + f + r.*aa - c;
    g = reshape(gg,T,2);


    return OrderedDict(:V => v , :r => r, :g => g, :adot => adot, :da => da, :astar => astar) 
end


# separately solve for the distribution => distribution conditional on optimal responses
function solve_FP(T,L,i_fix,a,da)

    # FOKKER-PLANCK EQUATION %
    LT = L';
    b = zeros(2T,1)

    #need to fix one value, otherwise matrix is singular
    i_fix = floor(Int64,T/2);
    b[i_fix]=.1;

    for j=1:T*2
	    LT[i_fix,j]=0;
    end

    LT[i_fix,i_fix]=1;

    #Solve linear system
    gg = LT\b;
    g_sum = gg'*ones(2T,1)*da;
    gg = gg./g_sum;
    g = reshape(gg,T,2);

    # aggregate capital
    KS = sum(g.*a*da);

    return OrderedDict(:KS => KS, :g => g)
end




