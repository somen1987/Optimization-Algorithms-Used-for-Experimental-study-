%  Sine Cosine Algorithm (SCA)  
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2011b(7.13)                                                                   
%                                                                                                     
%  Author and programmer: Seyedali Mirjalili                                                          
%                                                                                                     
%         e-Mail: ali.mirjalili@gmail.com                                                             
%                 seyedali.mirjalili@griffithuni.edu.au                                               
%                                                                                                     
%       Homepage: http://www.alimirjalili.com                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run SCA: [Best_score,Best_pos,cg_curve]=SCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________
function out=SCA(problem,params)
CostFunction=problem.CostFunction;

%fhd=str2func('cec13_func');

%load('extremum2013');


N=50;
dim=problem.nVar;


%target_err=1e-8;%when target error is greater thn 10^-8 then prog will terminate.
TEST_RUN=10;
fes=zeros(TEST_RUN,1);%initialize fes by zeros 
best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=5000;

Max_iteration=params.MaxIt;

lb=[problem.TxMin problem.TyMin problem.RotMin];
ub=[problem.TxMax problem.TyMax problem.RotMax];

for i=1:1:N
    X(i,:)=lb+(ub-lb).*rand(1,dim);
end
%Start_fun=input('Starting function#');
%End_fun=input('Ending function#');

disp('SCA is now tackling your problem')

%for func_num=Start_fun:1:End_fun

%tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
%filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

%load(filename);

%fprintf('Initial population is loaded\n');
%function [Destination_fitness,Destination_position,Convergence_curve]=SCA(N,Max_iteration,lb,ub,dim,fobj)

%display('SCA is optimizing your problem');

%Initialize the set of random solutions

% Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));
%fprintf('Calculate the fitness of the first set and find the best one\n');
% Calculate the fitness of the first set and find the best one
fitcount=0;

for i=1:size(X,1)
    Objective_values(1,i)=feval(CostFunction, X(i,:)');
    %fprintf('i=%d\n', i);
    fitcount=fitcount+1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    %All_objective_values(1,i)=Objective_values(1,i);
end

%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitnes
%err=abs(extremum(func_num) - Destination_fitness);
%fprintf('Initialization is completed\n');
while fitcount < FES %&& %err > target_err 
    
    % Eq. (3.4)
    a = 2;
%     Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Objective_values(1,i)= feval(CostFunction,X(i,:)');
        fitcount=fitcount+1;
        % Calculate the objective values
        
        
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    
    %Convergence_curve(t)=Destination_fitness;
    
    % Display the iteration and best optimum obtained so far
%     if mod(t,50)==0
%         display(['At iteration ', num2str(t), ' the optimum is ', num2str(Destination_fitness)]);
%     end
%     
    %Convergence_curve(t)=Destination_fitness;
    % Increase the iteration counter
    %t=t+1;
    
     %err=abs(extremum(func_num) - Destination_fitness);
end
 
   
    
   fprintf('global best=%f\n',Destination_fitness);
 
% best_it(test,1)=err;
 
 
 %fprintf('RUN %d is completed=================\n',test);
 
 fes(test,1)=fitcount;
end

%cpu_time=toc;
%T2=cpu_time/TEST_RUN

%best_run=min( best_it)
%worst_run=max( best_it )
%mean_run=mean( best_it)
%median_run=median(best_it);
%std_run=std( best_it)

%mean_fes=mean(fes)
%std_fes=std(fes)

%Succes_Ratio=(sum(best_it < 10e-8)/TEST_RUN)*100;
%filename='SCA.xlsx';
%sheet=sprintf('%dD',dim);
%if func_num>=26
    %column=sprintf('A%c','A'+func_num-26);
%else
 %    column=sprintf('%c','A'+func_num);
%end

%range=sprintf('%s2:%s51',column,column);

%xlswrite(filename,  best_it, sheet, range);
%xlswrite(filename,  best_it, sheet, range);

%range=sprintf('%s53',column);
%xlswrite(filename, best_run, sheet, range);

%range=sprintf('%s54',column);
%xlswrite(filename, worst_run, sheet, range);

%range=sprintf('%s55',column);
%xlswrite(filename, median_run, sheet, range);

%range=sprintf('%s56',column);
%xlswrite(filename, mean_run, sheet, range);

%range=sprintf('%s57',column);
%xlswrite(filename, std_run, sheet, range);

%range=sprintf('%s58',column);
%xlswrite(filename, mean_fes, sheet, range);

%range=sprintf('%s59',column);
%xlswrite(filename, std_fes, sheet, range);

%range=sprintf('%s61',column);
%xlswrite(filename, Succes_Ratio, sheet, range);
out.BestSol=Destination_position;
out.BestCost=Destination_fitness;
end%func

