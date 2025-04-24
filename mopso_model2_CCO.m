 %% MOPSO with Deb's Constraint Handling
clc; clear; close all;

%% Problem Definition
nVar = 6;                           % Number of Decision Variables
VarSize = [1 nVar];                % Size of Decision Variables Vector
VarMin = [0 0 0 0 0 0];            % Lower Bound
VarMax = [10000 10000 10000 10000 10000 10000];       % Upper Bound

CostFunction = @(x) EvaluateObjectives(x);
ConstraintFunction = @(x) EvaluateConstraints(x);

nObj = 2;

%% MOPSO Parameters
MaxIt = 100;           % Maximum Number of Iterations
nPop = 100;             % Population Size (Swarm Size)
nRep = 100;             % Repository Size

w = 0.5;               % Inertia Weight
c1 = 1.5;              % Personal Acceleration Coefficient
c2 = 2.0;              % Social Acceleration Coefficient

nGrid = 7;             % Number of Grids per Dimension
beta = 2;              % Leader Selection Pressure
gamma = 2;             % Deletion Selection Pressure

mu = 0.1;              % Mutation Rate

VelMax = 0.1 * (VarMax - VarMin);
VelMin = -VelMax;

%% Initialization
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Constraints = [];
empty_particle.Violation = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
empty_particle.Best.Violation = [];

particle = repmat(empty_particle, nPop, 1);

for i = 1:nPop
    particle(i).Position = VarMin + rand(VarSize) .* (VarMax - VarMin);
    particle(i).Velocity = zeros(VarSize);
    particle(i).Cost = CostFunction(particle(i).Position);
    [g, h] = ConstraintFunction(particle(i).Position);
    particle(i).Violation = ComputeViolation(g, h);
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.Violation = particle(i).Violation;
end

Repository = UpdateRepository(particle, nRep);

%% MOPSO Main Loop
for it = 1:MaxIt
    for i = 1:nPop
        leader = SelectLeader(Repository, beta);
        particle(i).Velocity = ...
            w * particle(i).Velocity ...
            + c1 * rand(VarSize) .* (particle(i).Best.Position - particle(i).Position) ...
            + c2 * rand(VarSize) .* (leader.Position - particle(i).Position);

        particle(i).Velocity = max(particle(i).Velocity, VelMin);
        particle(i).Velocity = min(particle(i).Velocity, VelMax);

        particle(i).Position = particle(i).Position + particle(i).Velocity;
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);

        particle(i).Cost = CostFunction(particle(i).Position);
        [g, h] = ConstraintFunction(particle(i).Position);
        particle(i).Violation = ComputeViolation(g, h);

        % Apply Mutation
        if rand < mu
            particle(i).Position = Mutate(particle(i).Position, mu, VarMin, VarMax);
        end

        if Dominates(particle(i), particle(i).Best)
            particle(i).Best = particle(i);
        elseif ~Dominates(particle(i).Best, particle(i))
            if rand < 0.5
                particle(i).Best = particle(i);
            end
        end
    end

    Repository = UpdateRepository(particle, nRep);
    disp(['Iteration=' num2str(it)]);
end

%% Results
PlotCosts(Repository);



%% Objective Function 
function z = EvaluateObjectives(x)
    % Decision variables
    k = round(x(1));
    l = round(x(2));
    m = round(x(3));
    n = round(x(4));
    q = x(5);
    T = x(6);

  % Parameters
    d=50000; s_m=400; d_m=120; gamma_m=24000; gamma_r=30000; c_m=150; s_mr=25;
    sigma=140000; lambda=160000; hm1=3; l1=0.074; l2=0.074;
    E_delta=0.8; hm2=4; A_r=300; w_r=150; s_r_dash=15000; I_r=175200;
    h_r1=2; h_r2=3; tau=1; h_c=5; l3=0.074; E_L=0; alpha=0.5;
    c_r=5; beta=0.8; x_param=8500;  a_re=3*10^(-7); b_re=0.0012; c_re=1.4;
    a_me=3*10^(-7); b_me=0.0012; c_me=1.4; L_mr=100; L_rc=100; L_cm=200;
    E=36.75; v_e=2.6; w_e=10; d_e=100;

  % Carbon policy parameters
    C_max = 1000;  % Example value, update as needed
    c_t = 5;       % Example penalty cost, update as needed

    %for the first cycle
    trk1=0; tpk1=0;u_1=0;


   
    % Precompute summations
    A1 = 0; A2 = 0; A3 = 0; A4 = 0;

    for j = 1:k
        A1 = A1 + ((1+trk1+j*T)^(2+l1)-(1+trk1+(j-1)*T)^(2+l1));
        A2 = A2 + ((1+trk1+(j-1)*T)^(1+l1));
    end
    for j = 1:k-1
        A3 = A3 + ((1+trk1+j*T)^(1+l1)-(1+trk1)^(1+l1));
    end
    for j = l+1:k
        A4 = A4 + ((1+tpk1+(j-k+1)*T)^(2+l2)-(1+tpk1+(j-k)*T)^(2+l2));
    end

    % Objective Function 1 (base version without carbon penalty)
    f_1 = (d/(l*q+(m-l)*(1-E_delta)*q)) * ...
        (s_m+(m-l)*q*d_m+(m-l)*T*gamma_m+l*T*gamma_r+(m-l)*E_delta*q*c_m+ ...
        m*q*s_mr+(sigma/((1+l1)*(2+l2))*A1-(sigma*T)/(1+l1)*A2+(sigma*T)/(1+l1)*A3+ ...
        (l-k)*(sigma*T)/(1+l1)*((1+trk1+k*T)^(1+l1)-(1+trk1)^(1+l1))- ...
        0.5*l*(l-1)*q*T)*hm1 + ...
        (lambda/((1+l2)*(2+l2))*((1+tpk1+(l-k+1)*T)^(2+l2)-(1+tpk1)^(2+l2))- ...
        (1/(1+l2))*(l-k+1)*lambda*T*(1+tpk1)^(1+l2)+ ...
        lambda/((1+l2)*(2+l2))*A4- ...
        0.5*(n^2-n-l^2+l)*(1+tpk1)^(1+l2)+ ...
        T*(m-n)*lambda/(1+l2)*((1+tpk1+(n-k)*T)^(1+l2)-(1+tpk1)^(1+l2))- ...
        0.5*(m-l-1)*q*T)*hm2 + ...
        A_r + m*q*w_r+(m-l)*q*s_r_dash/I_r+ ...
        0.5*m*(m+1)*E_delta*q*T*h_r1 + ...
        (0.5*l^2*(q-d*T)*T + 0.5*l*q*T + l*(m-l)*T*(q-d*T) + ...
        0.5*(m-l+1)*(m-l)*T*((1-E_delta)*q-d*T)+(m-l)*q*T- ...
        (m-l)*(T-tau)*(E_delta*q+d*tau)- ...
        0.5*(m-l)*d*tau^2-0.5*(m-l)*d*(T-tau)^2+ ...
        0.5*d*((m+l)^2*(q-d*T)-2*(m+1)*(m-l+1)*E_delta*q*(q-d*T)+(m-l)^2*E_delta^2*q^2)) * h_r2 + ...
        alpha*(l*q+(m-l)*(1-E_delta)*q)*c_r + ...
        1/d*(l^2*q^2+2*l*(m-l)*(1-E_delta)*q+(m-l)^2*(1-E_delta)^2*q^2)*h_c + ...
        ((1-beta)*E_L*(1/((1-l3)*x_param))*((u_1+E_L)^(1-l3)-u_1^(1-l3)) + ...
        beta*E_L*k*T - ...
        (sigma/((1+l1)*(2+l1)))*((1+trk1+k*T)^(2+l1)-(1+trk1)^(2+l1)) - ...
        (sigma*k*T)/(1+l1)*(1+trk1)^(1+l1)) * h_c);

    % Objective Function 2
    f2 = (d/(l*q+(m-l)*(1-E_delta)*q)) * ...
        (a_re*(sigma^2/(1+2*l1))*((1+trk1+k*T)^(1+2*l1)-(1+trk1)^(1+2*l1)) - ...
        b_re*sigma*(1+l1)*((1+trk1+k*T)^(1+l1)-(1+trk1)^(1+l1)) + ...
        c_re*k*T + ...
        a_me*(lambda^2/(1+2*l2))*((1+tpk1+(n-k)*T)^(1+2*l2)-(1+tpk1)^(1+2*l2)) - ...
        b_me*lambda*(1+l2)*((1+tpk1+(n-k)*T)^(1+l2)-(1+tpk1)^(1+l2)) + ...
        c_me*(n-k)*T + ...
        (2*m*L_mr+2*L_rc+2*L_cm)*E*v_e + ...
        (0.5*m*(m+1)*E_delta*q*T + ...
        1/d*(l^2*q^2+2*l*(m-l)*(1-E_delta)*q+(m-l)^2*(1-E_delta)^2*q^2))*w_e + ...
        (1-beta)*E_L*d_e);

    % Apply carbon cap and offset policy
    if f2 > C_max
        f1 = f_1 + c_t * (f2 - C_max);
    else
        f1 = f_1;
    end

    z = [f1, f2];
end




%% Constraints
function [g, h] = EvaluateConstraints(x)
    % Defining the constraints
    k = round(x(1));
    l = round(x(2));
    m = round(x(3));
    n = round(x(4));
    q = x(5);
    T = x(6);
    
   % Parameters
    d=50000; sigma=140000; lambda=160000; trk1=0; tpk1=0;  l1=0.074; l2=0.074;
    E_delta=0.8;
    E_L=0; 
   beta=0.8; 
   

    % Constraints
    g(1) = l * q - (sigma / (1 + l1)) * ((1 + trk1 + k * T)^(1 + l1) - (1 + trk1)^(1 + l1)); % g1 <= 0
    g(2) = (m - l) * q - (lambda / (1 + l2)) * ((1 + tpk1 + (n - k) * T)^(1 + l1) - (1 + tpk1)^(1 + l2)); % g2 <= 0
    g(3) = q - (sigma / (1 + l1)) * ((1 + trk1 + T)^(1 + l1) - (1 + trk1)^(1 + l1)); % g3 <= 0
    g(4) = q - (lambda / (1 + l2)) * ((1 + tpk1 + T)^(1 + l1) - (1 + tpk1)^(1 + l2)); % g4 <= 0
    g(5) = d * T - q - E_delta * q; % g5 <= 0
    h(1) = (((1 + l1) * beta / sigma) * E_L + (1 + trk1)^(1 + l1)) - (sigma * T) / (1 + 2 * l1); % h1 = 0
    
end

%% Constraint Violation
function v = ComputeViolation(g, h)
    tol = 1e-4;
    g_violation = sum(max(g, 0));
    h_violation = sum(abs(h) > tol .* abs(h));
    v = g_violation + h_violation;
end

%% Domination Check
function flag = Dominates(x, y)
    if x.Violation < y.Violation
        flag = true;
    elseif x.Violation > y.Violation
        flag = false;
    else
        flag = all(x.Cost <= y.Cost) && any(x.Cost < y.Cost);
    end
end

%% Repository Update
function rep = UpdateRepository(pop, nRep)
    rep = pop;
    to_remove = false(size(rep));

    for i = 1:numel(rep)
        for j = 1:numel(rep)
            if i ~= j && Dominates(rep(j), rep(i))
                to_remove(i) = true;
                break;
            end
        end
    end

    rep = rep(~to_remove);

    if numel(rep) > nRep
        Costs = vertcat(rep.Cost);
        D = pdist2(Costs, Costs);
        CrowdingDistance = sum(D, 2);
        [~, sortIdx] = sort(CrowdingDistance, 'descend');
        rep = rep(sortIdx(1:nRep));
    end
end

%% Leader Selection
function leader = SelectLeader(rep, beta)
    n = numel(rep);
    P = rand(1, n).^beta;
    [~, idx] = min(P);
    leader = rep(idx);
end

%% Mutation
function y = Mutate(x, mu, VarMin, VarMax)
    nVar = numel(x);
    nmu = ceil(mu * nVar);
    j = randperm(nVar, nmu);
    sigma = 0.1 * (VarMax - VarMin);
    y = x;
    y(j) = x(j) + sigma(j) .* randn(size(j));
    y = max(y, VarMin);
    y = min(y, VarMax);
end


