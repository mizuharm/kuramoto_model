%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kuramoto model with inertia dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Simulate dynamics for single set of parameters
%
 

function [t_vec,data_matrix] = solve_kuramoto_ode(uint,a,w,k,n,G,endTime,dt,model,alpha)

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve ODE
    %%%%%%%%%%%%%%%%%%%%%%%%%%
   
        %opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
        [t,u]=ode45(@(t,y) kuramoto_inertia(y,a,w,k,n,G,model,alpha),[0,endTime],[uint(1:n),uint(n+1:2*n)]);
       
        %Interpolate time points for even movie 

        t_vec = 0:dt:endTime;
        u = interp1(t,u,t_vec);
        u(:,1:n) = mod(u(:,1:n),2*pi);

        data_matrix = u; 
end

function output = kuramoto_inertia(u,a,w,k,n,G,model,alpha)
    %Function : defines the Kuramoto model with inertia
    %           to be used in ode45
 
    %Inputs: n - number of oscillators 
    %        u - vector of length 2n; oscillator phases and velocities
    %        a - drag coefficient
    %        k - coupling strength
    %        G - adjacency matrix
    
    pos = u(1:n);
    vel = u(n+1:end);
    u_mat = repmat(pos,1,n);
    F = (k/n)*sum(G.*sin(u_mat' - u_mat-alpha),2);
    if model==2 %second order model
        output = [vel; -a*vel + w + F];
    elseif model==1 %first order model
        output = [w+F;-vel+w+F];
    end
end