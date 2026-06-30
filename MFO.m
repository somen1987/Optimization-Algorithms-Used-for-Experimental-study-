%______________________________________________________________________________________________
%  Moth-Flame Optimization Algorithm (MFO) toolbox                                                            
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
%  S. Mirjalili, Moth-Flame Optimization Algorithm: A Novel Nature-inspired Heuristic Paradigm, 
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.07.006
%_______________________________________________________________________________________________
% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run MFO: [Best_score,Best_pos,cg_curve]=MFO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________

function out=MFO(problem,params)
global fixedimg transformedimg
CostFunction=problem.CostFunction;

%fhd=str2func('cec13_func');

%load('extremum2013');


N=50;
%dim=input('Dimension:');
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
%Initialize the positions of moths
%Moth_pos=initialization(N,dim,ub,lb);

%for func_num=Start_fun:1:End_fun

tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
%filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

%load(filename);

Moth_pos=X;

Convergence_curve=zeros(1,Max_iteration);

Iteration=1;
fitcount=0;
% Main loop
%while Iteration<Max_iteration+1
%err=1;
while fitcount < FES %&& err > target_err 

    
    % Number of flames Eq. (3.14) in the paper
    Flame_no=round(N-Iteration*((N-1)/Max_iteration));
    
    for i=1:size(Moth_pos,1)
        
        % Check if moths go out of the search spaceand bring it back
        Flag4ub=Moth_pos(i,:)>ub;
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        % Calculate the fitness of moths
        Moth_fitness(1,i)=feval(CostFunction, Moth_pos(i,:)');
        fitcount=fitcount+1;
        All_fitness(1,i)=Moth_fitness(1,i);
        
    end
       
    if Iteration==1
        % Sort the first population of moths
        [fitness_sorted I]=sort(Moth_fitness);
        sorted_population=Moth_pos(I,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    else
        
        % Sort the moths
        double_population=[previous_population;best_flames];
        double_fitness=[previous_fitness best_flame_fitness];
        
        [double_fitness_sorted I]=sort(double_fitness);
        double_sorted_population=double_population(I,:);
        
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    end
    
    % Update the position best flame obtained so far
    Best_flame_score=fitness_sorted(1);
    Best_flame_pos=sorted_population(1,:);
      
    previous_population=Moth_pos;
    previous_fitness=Moth_fitness;
    
    % a linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a=-1+Iteration*((-1)/Max_iteration);
    
    for i=1:size(Moth_pos,1)
        
        for j=1:size(Moth_pos,2)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                
                % D in Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
            end
            
            if i>Flame_no % Upaate the position of the moth with respct to one flame
                
                % Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
            end
            
        end
        
    end
    
    %Convergence_curve(Iteration)=Best_flame_score;
    
%      if Iteration>2
%         line([Iteration-1 Iteration], [Convergence_curve(Iteration-1) Convergence_curve(Iteration)],'Color','b')
%         xlabel('Iteration');
%         ylabel('Best score obtained so far');        
%         drawnow
%     end
 
    
%     set(handles.itertext,'String', ['The current iteration is ', num2str(Iteration)])
%     set(handles.optimumtext,'String', ['The current optimal value is ', num2str(Best_flame_score)])
%     if value==1
%         hold on
%         scatter(Iteration*ones(1,N),All_fitness,'.','k')
%     end
    
    Iteration=Iteration+1; 
    
    %err=abs(extremum(func_num) - Best_flame_score);
    
end

fprintf('global best=%f\n',Best_flame_score);
 
 %best_it(test,1)=err;
 
 
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
%filename='MFO.xlsx';
%sheet=sprintf('%dD',dim);
%if func_num>=26
    %column=sprintf('A%c','A'+func_num-26);
%else
     %column=sprintf('%c','A'+func_num);
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
out.BestSol=Best_flame_pos;
out.BestCost=Best_flame_score;
image_transformation(out.BestSol);

%imtool(transformedimg)
out.BestPcc=pcc(fixedimg,transformedimg);
end%func



