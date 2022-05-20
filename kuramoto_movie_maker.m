%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kuramoto model with inertia dynamics movie maker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Simulates Kuramoto model with inertia for single set of parameters
%
%    Model: u_i'' +a u_i' = w_i + K/n * sum(G_{ij} sin(u_j-u_i-alpha)
%           i = {1,...,n}
%    
%    Code contents:
%       -User chooses system parameters
%       -System solves Kuramoto model from random IC and plots data
%       -Outputs MP4 movie file of plots
%
%
%        Inputs:
%         n = number of oscillators
%        model = order of Kuramoto Model
             %model = 1; %First order
             %model = 2; %Second order
%        a = damping coefficient for second order model
%        k = coupling strength
%        alpha = phase lag
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kuramoto_movie_maker(n,model,a,k,alpha)

        tic %track total runtime 
        
        %add auxiliary files
        addpath([pwd, '\movie_maker_aux_files'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        %Set total runtime
            endTime = 200;

        %Timestep in movie (interpolates data from ode45)
            dt = .3;
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Adjacency Matrix
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     disp('What is the graph coupling?')
     
     adj_matrix_choice = input('Options: "Nearest Neighbor", "Erdos-Renyi", "All-to-All"\n','s');
     
     matrix_chosen = 0;
     
     while matrix_chosen == 0
         if adj_matrix_choice == "Nearest Neighbor"
            %NN coupling range
            r = input('What is the coupling range for Nearest Neighbor? (0<r<.5)\n');
            G = nn_graph(n,r);
            matrix_chosen = 1;
         elseif adj_matrix_choice == "Erdos-Renyi"
            %ER coupling probability
            p = input('What is the connection probability for Erdos-Renyi? (0<p<1)\n');      
            G = er_graph(n,p);
            matrix_chosen = 1;
         elseif adj_matrix_choice == "All-to-All"
            G=ones(n,n)-eye(n,n);
            matrix_chosen = 1;
         else
             %Non-valid choice
             disp('Option not valid; please input another coupling choice.')
             adj_matrix_choice = input('Options: "Nearest Neighbor", "Erdos-Renyi", "All-to-All"\n','s');
         end
     end
 

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Distribution Choice
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     disp('What is the intrinsic frequency distribution?')
     
     distribution_type = input('Options: "Unimodal Gaussian", "Bimodal Gaussian"\n','s');
       
     distribution_chosen = 0;
     
     while distribution_chosen == 0
         if distribution_type == "Unimodal Gaussian"
             mu = input('What is the mean of the distribution?\n');
             sig = input('What is the standard deviation of the distribution? \n');
             w = normrnd(mu,sig,[n,1]);
             distribution_chosen = 1;
         elseif distribution_type == "Bimodal Gaussian"
            mu1 = input('Where is the first peak of the distribution?\n');
            sig1 = input('What is the standard deviation of the first peak?\n');
            mu2 = input('Where is the second peak of the distribution?\n');
            sig2 =  input('What is the standard deviation of the second peak?\n');
   
            %Assume equal mass of peaks; can manually change
            size1 = floor(.5*n);
            size2 =  n-size1;

            w1= normrnd(mu1,sig1,[size1,1]);
            w2 = normrnd(mu2,sig2,[size2,1]);
            w=[w1;w2];
            distribution_chosen = 1;
         else 
            %Non-valid choice
            disp('Option not valid; please input another distribution choice.')   
            distribution_type = input('Options: "Unimodal Gaussian", "Bimodal Gaussian"\n','s');
         end
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initial conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Phases uniformly around circle
        u_phase_int= 2*pi*rand(n,1);
        
        %Velocities chosen as intrinsic frequencies
        u_prime_int = w;

        uint = [u_phase_int;u_prime_int];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve the Kuramoto equation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [t,u] = solve_kuramoto_ode(uint,a,w,k,n,G,endTime,dt,model,alpha);
 
        %Calculate Order Parameter
            ord_param = mean(exp(1i*u(:,1:n)),2);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot results and save to movie
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         fig=figure();

        %Figure size options
 
            set(gcf, 'Position',  [100 206 1100 430])
 
        for i =1:length(t)
            clf

            subplot(1,2,1)
                plot(1:n,u(i,1:n),'.','MarkerSize',16)
                a=gca();
                axis([1 n 0 2*pi])
                a.FontSize = 18;
                a.XTick = [1 n/2 n];
                a.YTick = [0, pi, 2*pi];
                a.YTickLabel = {'0','\pi','2\pi'}; 
                xlabel('$i$','interpreter','latex')
                ylabel('$u_i$','interpreter','latex')
 
            subplot(1,2,2)
                X = 0:.01:2*pi; %plot unit circle
                hold on 
                plot(cos(X),sin(X))
                plot(cos(u(i,1:n)),sin(u(i,1:n)),'.','MarkerSize',16)
                plot([0,real(ord_param(i))],[0,imag(ord_param(i))],'-','LineWidth',3)
                hold off
                a=gca();
                axis([-1 1 -1 1])
                a.FontSize = 18;
                a.XTick = [-1 0 1];
                a.YTick = [-1, 0, 1];
            
            F(i) = getframe(fig);
        end

        v=VideoWriter('kuramoto.mp4','MPEG-4');
        open(v);
        writeVideo(v,F);
        close(v); 

    toc
    
end

