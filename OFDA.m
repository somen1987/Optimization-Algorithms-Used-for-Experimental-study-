% Opposition Flow Direction Algorithm (OFDA) source Code Version 1.0
%
% Developed in MATLAB R2018b
%
% Author and programmer: 
% Dr Manoj Kumar Naik
% Faculty of Engineering and Technology, Siksha O Anusandhan, Bhubaneswar, Odisha-751030, India 
% e-mail:       naik.manoj.kumar@gmail.com
% ORCID:        https://orcid.org/0000-0002-8077-1811
% SCOPUS:       https://www.scopus.com/authid/detail.uri?authorId=35753522900
% WOS:          https://www.webofscience.com/wos/author/record/O-2982-2017
% G-Scholar:    https://scholar.google.co.in/citations?user=tX-8Xw0AAAAJ&hl=en 
% Researchgate: https://www.researchgate.net/profile/Manoj_Naik9
% DBLP:         https://dblp.uni-trier.de/pers/k/Kumar:Naik_Manoj
%_____________________________________________________________________________________________________           
% Please cite to the main paper:
% R. Panda, M. Swain, M. K. Naik, S. Agrawal, and A. Abraham, 
% �A Novel Practical Decisive Row-class Entropy-based Technique for Multilevel Threshold Selection 
% Using Opposition Flow Directional Algorithm,� IEEE Access, p. 1, 2022, doi: 10.1109/ACCESS.2022.3215082.
%
% This program using the framework of FDA by SEYEDALI MIRJALILI
% https://seyedalimirjalili.com/projects
%_____________________________________________________________________________________________________
function out=OFDA(problem,params)
% Initialize the positions of flows
global fixedimg transformedimg;
alpha=30;
maxiter=params.MaxIt;
lb=[problem.TxMin problem.TyMin problem.RotMin];
ub=[problem.TxMax problem.TyMax problem.RotMax];
dim=3;
fobj=problem.CostFunction;
beta=1; %Nighbourhood points
flow_x=initialization(alpha,dim,ub,lb);
neighbor_x=zeros(beta,dim);
newflow_x=inf(size(flow_x));
newfitness_flow=inf(size(flow_x,1));
ConvergenceCurve=zeros(1,maxiter);
fitness_flow=inf.*ones(alpha,1);
fitness_neighbor=inf.*ones(beta,1);
%% calculate fitness function of each flow
for i=1:alpha
    fitness_flow(i,:)=fobj(flow_x(i,:));%fitness of each flow
end
%% sort results and select the best results
[~,indx]=sort(fitness_flow);
flow_x=flow_x(indx,:);
fitness_flow=fitness_flow(indx);
Best_fitness=fitness_flow(1);
BestX=flow_x(1,:);
%% Initialize velocity of flows
Vmax=0.1*(ub-lb);
Vmin=-0.1*(ub-lb);
%% Main loop
for iter=1:maxiter
    % Update W
    %ofitness_flow=fitness_flow;
    W=(((1-1*iter/maxiter+eps)^(2*randn)).*(rand(1,dim).*iter/maxiter).*rand(1,dim));
    % Update the Position of each flow
    for i=1:alpha
        % Produced the Position of neighborhoods around each flow
        for j=1:beta
            Xrand=lb+rand(1,dim).*(ub-lb);
            delta=W.*(rand*Xrand-rand*flow_x(i,:)).*norm(BestX-flow_x(i,:));
            neighbor_x(j,:)=flow_x(i,:)+randn(1,dim).*delta;
            neighbor_x(j,:)=max(neighbor_x(j,:),lb);
            neighbor_x(j,:)=min(neighbor_x(j,:),ub);
            fitness_neighbor(j)=fobj(neighbor_x(j,:));
        end
        % Sort position of neighborhoods
          [~,indx]=sort(fitness_neighbor);
          % Update position, fitness and velocity of current flow if the fitness of best neighborhood is
          % less than of current flow
          if fitness_neighbor(indx(1))<fitness_flow(i)
              % Calculate slope to neighborhood
              Sf=(fitness_neighbor(indx(1))-fitness_flow(i))./sqrt(norm(neighbor_x(indx(1),:)-flow_x(i,:)));%calculating slope
              % Update velocity of each flow
              V=randn.*(Sf);
              if V<Vmin
                  V=-Vmin;
              elseif V>Vmax
                  V=-Vmax;
              end
              %Flow moves to best neighborhood
              newflow_x(i,:)=flow_x(i,:)+V.*(neighbor_x(indx(1),:)-flow_x(i,:))./sqrt(norm(neighbor_x(indx(1),:)-flow_x(i,:)));
          else
              %Generate integer random number (r)
              r=randi([1 alpha]);
              % Flow moves to r th flow if the fitness of r th flow is less
              % than current flow
             if fitness_flow(r)<=fitness_flow(i)
                 newflow_x(i,:)=flow_x(i,:)+randn(1,dim).*(flow_x(r,:)-flow_x(i,:));
              else
                 newflow_x(i,:)=flow_x(i,:)+randn*(BestX-flow_x(i,:));
             end
          end
          % Return back the flows that go beyond the boundaries of the search space
              newflow_x(i,:)=max(newflow_x(i,:),lb);
              newflow_x(i,:)=min(newflow_x(i,:),ub);
         % Calculate fitness function of new flow 
              newfitness_flow(i)=fobj(newflow_x(i,:));
         % Update current flow     
          if newfitness_flow(i)<fitness_flow(i)
              flow_x(i,:)=newflow_x(i,:);
              fitness_flow(i)=newfitness_flow(i);
          end
         % Update  best flow 
         if fitness_flow(i)<Best_fitness
             BestX=flow_x(i,:);
             Best_fitness=fitness_flow(i);
         end
    end
     %% Oppsotion based learning
     for i=1:alpha
            D(i,:)=min(flow_x(i,:))+max(flow_x(i,:))-flow_x(i,:);
            Flag4ub=D(i,:)>ub;
            Flag4lb=D(i,:)<lb;
            D(i,:)=(D(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            op_fitness=fobj(D(i,:));
            if op_fitness<fitness_flow(i)
                flow_x(i,:)=D(i,:);
                fitness_flow(i)=op_fitness;
            end

     end
     
     %%
    ConvergenceCurve(iter)=Best_fitness;
    MeanPos = mean(flow_x);
div = 0;
for k = 1:size(flow_x,1)
    div = div + norm(flow_x(k,:) - MeanPos);
end
DiversityCurve(iter) = div / size(flow_x,1);
    disp(['MaxIter= ' ,num2str(iter), 'BestFit= ', num2str(Best_fitness)])%disply results
end

out.BestSol=BestX;
out.BestCost=Best_fitness;
out.ConvergenceCurve=ConvergenceCurve;
out.DiversityCurve=DiversityCurve;
image_transformation(out.BestSol);

%imtool(transformedimg)
out.BestPcc=pcc(fixedimg,transformedimg);
out.Bestssim=ssim(fixedimg,transformedimg);
end