function out=PSO(problem,params)
global fixedimg movingimg transformedimg
%% Problem Definition
objfun = problem.CostFunction;  % Example function (Sphere)

dim = problem.nVar;
%lb  = -10 * ones(1, dim);
%ub  =  10 * ones(1, dim);
lb=[problem.TxMin problem.TyMin problem.RotMin];% Lower Bound of Decision Variables
ub=[problem.TxMax problem.TyMax problem.RotMax]; % Upper Bound of Decision Variables

%% PSO Parameters
pop_size = 50;
max_iter = 100;

w  = 0.7;    % inertia weight
c1 = 1.5;    % cognitive coefficient
c2 = 1.5;    % social coefficient

%% Initialization
pos = rand(pop_size, dim) .* (ub - lb) + lb;
vel = zeros(pop_size, dim);

fitness = zeros(pop_size,1);

for i = 1:pop_size
    fitness(i) = objfun(pos(i,:));
end

pbest = pos;
pbest_fit = fitness;

[gbest_fit, idx] = min(fitness);
gbest = pos(idx,:);

convergence = zeros(1, max_iter);
diversity   = zeros(1, max_iter);

%% Main Loop
for iter = 1:max_iter
    
    for i = 1:pop_size
        
        % Velocity update
        vel(i,:) = w * vel(i,:) ...
                 + c1 * rand(1,dim) .* (pbest(i,:) - pos(i,:)) ...
                 + c2 * rand(1,dim) .* (gbest - pos(i,:));
        
        % Position update
        pos(i,:) = pos(i,:) + vel(i,:);
        
        % Boundary control
        pos(i,:) = max(min(pos(i,:), ub), lb);
        
        % Evaluate
        fitness(i) = objfun(pos(i,:));
        
        % Update personal best
        if fitness(i) < pbest_fit(i)
            pbest(i,:) = pos(i,:);
            pbest_fit(i) = fitness(i);
        end
        
        % Update global best
        if fitness(i) < gbest_fit
            gbest = pos(i,:);
            gbest_fit = fitness(i);
        end
    end
    
    convergence(iter) = gbest_fit;
    
    %% 🔥 Population Diversity
    mean_pop = mean(pos, 1);
    div = 0;
    
    for i = 1:pop_size
        for j = 1:dim
            div = div + abs(pos(i,j) - mean_pop(j)) / (ub(j) - lb(j));
        end
    end
    
    diversity(iter) = div / (pop_size * dim);
    
    fprintf('Iter %d | Best = %f | Diversity = %f\n', ...
        iter, gbest_fit, diversity(iter));
end

%% 📊 Plot Results
% figure;
% 
% subplot(2,1,1);
% plot(convergence, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Best Fitness');
% title('PSO Convergence Curve');
% grid on;
% 
% subplot(2,1,2);
% plot(diversity, 'r', 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Population Diversity');
% title('PSO Diversity Curve');
% grid on;

%% Results
disp('Best Solution:');
disp(gbest);

disp('Best Fitness:');
disp(gbest_fit);

out.BestSol=gbest;
out.BestCost=gbest_fit;
out.convergence=convergence;
out.diversity=diversity;

image_transformation(out.BestSol);

%imtool(transformedimg)
out.BestPcc=pcc(fixedimg,transformedimg);
image_transformation(out.BestSol);

%imtool(transformedimg)
out.Bestssim=ssim(fixedimg,transformedimg);
image_transformation(out.BestSol);
end