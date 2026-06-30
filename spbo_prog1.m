function out=spbo_prog1(problem,params)
CostFunction =problem.CostFunction;
%fhd=str2func('cec13_func');

%load('extremum2013');


%dim=input('Dimension:');


variable= problem.nVar;
student = 30;
N=student;

%target_err=1e-8;%when target error is greater thn 10^-8 then prog will terminate.
TEST_RUN=1;
fes=zeros(TEST_RUN,1);%initialize fes by zeros 
best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=5000;

%Max_iteration=round(FES/N);

Max_iteration=params.MaxIt;



mini=ones(1,variable).*(-100); %where mini is the lower bound of variable n
maxi=ones(1,variable).*100; %where maxi is the upper bound of variable n





%function [Best_fitness,Best_student,Convergence_curve]=SPBO(student,Max_iteration,maxi,mini,variable,fobj)

display('SPBO is optimizing your problem');




%Start_fun=input('Starting function#');
%End_fun=input('Ending function#');


%for func_num=Start_fun:1:End_fun

tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
%filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

%load(filename);

%fprintf('Initial population is loaded\n');



%Initialize the set of random solutions

for i=1:1:student
    X(i,:)=mini+(maxi-mini).*rand(1,variable);
end



sol=zeros(1,variable);
ans=inf;

%Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(size(X,1),1);
fitcount=0;
% Calculate the fitness of the first set and find the Best_fitness one
for i=1:student
    Objective_values(i,1)=CostFunction(X(i,:));
    fitcount=fitcount+1;
        
end
% display (sol);

% display (ans);




 %-----Opposition-based Learning
dynXmin=min(X);
dynXmax=max(X);

centroid=mean(X);

beta0=180;

for i=1:1:student
    
    U=X(i,:)-centroid;
    
    V=sqrt((X(i,:)-dynXmin).*(dynXmax-X(i,:)));
    
    Opp_X(i,:)=centroid+ U.*cos(beta0)-V.*sin(beta0);
    
    
    Opp_Fitness(i)=CostFunction(X(i,:));
    fitcount=fitcount+1;
end
        
        AllFitness=[Objective_values' Opp_Fitness];
        All_X=[X; Opp_X];
   
        [sorted_fitness index]=sort(AllFitness);
        X=All_X(index(1:student),:);
        Objective_values=sorted_fitness(1:student)';
        
        Best_fitness=sorted_fitness(1,1);
        Best_student=X(1,:);

        
        sol=X;
        ans=Objective_values;
        

   %err=abs(extremum(func_num) - Best_fitness);
Pgj=0.3;   
%Main loop   
%for t=1:1:Max_iteration
while fitcount < FES %&& err > target_err 
    
    
    if rand < Pgj
    
   for do=1:1:variable
        
%     sum1=zeros(1,variable);
%     for gw=1:1:variable
%     for fi=1:1:student
%         sum1(1,gw)=sum1(1,gw)+sol(fi,gw);
%     end;
%     pop_mean(1,gw)=sum1(1,gw)/student;
%     end;
   pop_mean=mean(sol);
  
     par=sol;
     par1=sol;
    
  
    
     check=rand(student,1);
    mid=rand(student,1);
    for dw=1:1:student
       % Best Student
        if Best_fitness==ans(dw,1)
            
             jg=ans(randperm(numel(ans),1));
         
         for oi=1:1:student
             if jg==ans(oi,1)
                 lk=oi;
             end;
         end;
     
            par1(dw,do)=par(dw,do)+(((-1)^(round(1+rand)))*rand*(par(dw,do)-par(lk,do)));       % Equation (1)
           
        else if check(dw,1)<mid(dw,1)
         % Good Student
                rta=rand;
                if rta>rand
                    par1(dw,do)=Best_student(1,do)+(rand*(Best_student(1,do)-par(dw,do)));      % Equation (2a)
                else
                
                par1(dw,do)=par(dw,do)+(rand*(Best_student(1,do)-par(dw,do)))+((rand*(par(dw,do)-pop_mean(1,do))));         % Equation (2b)
                end;
            else
                an=rand;
          % Average Student
                if rand>an
                    
                    par1(dw,do)=par(dw,do)+(rand*(pop_mean(1,do)-par(dw,do)));      % Equation (3)
                  
                else
           % Students who improves randomly
                        par1(dw,do)=mini(do)+(rand*(maxi(do)-mini(do)));                    % Equation (4)
                   
                end;
            end;
        end;
    end;

   
    % Boundary checking of the improvement of the students
    for z=1:1:student
       
            if par1(z,do)>maxi(do)
                par1(z,do)=maxi(do);
            else if par1(z,do)<mini(do)
                    par1(z,do)=mini(do);
                end;
            
        end;
    end;
    
    X=par1;
    
   for i=1:1:size(X,1)
        % Calculate the objective values
        Objective_values(i,1)=CostFunction(X(i,:));
        fitcount=fitcount+1;
   end;
        
        
       fun1=Objective_values;

        % Update the solution if there is a better solution
        for vt=1:1:student
        if ans(vt,1)>fun1(vt,1)
            ans(vt,1)=fun1(vt,1);
            sol(vt,:)=par1(vt,:);
        end;
    end;
       
     Best_fitness1=min(ans);
     for fo=1:1:student
             if Best_fitness1==ans(fo,1)
                 Best_student1=sol(fo,:);
             end;
         end;
         
         % Update the best student
          if Best_fitness>Best_fitness1
         Best_fitness=Best_fitness1;
         Best_student=Best_student1;
     end;
      end;
    else
 %-----Opposition-based Learning
dynXmin=min(X);
dynXmax=max(X);

centroid=mean(X);

beta0=180;

for i=1:1:student
    
    U=X(i,:)-centroid;
    
    V=sqrt((X(i,:)-dynXmin).*(dynXmax-X(i,:)));
    
    Opp_X(i,:)=centroid+ U.*cos(beta0)-V.*sin(beta0);
    
    
    Opp_Fitness(i)=CostFunction(X(i,:));
    fitcount=fitcount+1;
end
        
        AllFitness=[Objective_values' ;Opp_Fitness];
        All_X=[X; Opp_X];
   
        [sorted_fitness index]=sort(AllFitness);
        X=All_X(index(1:student),:);
        Objective_values=sorted_fitness(1:student)';
        
        Best_fitness1=sorted_fitness(1,1);
        Best_student1=X(1,:);
        
        if Best_fitness>Best_fitness1
         Best_fitness=Best_fitness1;
         Best_student=Best_student1;
        end
        
        
end
      
   % Display the iteration and Best_fitness optimum obtained so far 
%     Convergence_curve(t)=Best_fitness;
%     display (t);
%     display (Best_fitness);
    
    %err=abs(extremum(func_num) - Best_fitness);
    
end


fprintf('global best=%f \n',Best_fitness);
 
 %best_it(test,1)=err;
 
 
 fprintf('RUN %d is completed=================\n',test);
 
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

%display(Succes_Ratio)

%filename='SPBO_ROBL.xlsx';
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


%end%func
out.BestSol=Best_student;
out.BestCost=Best_fitness;
SPBO_ROBL_T2pre_T2post_best_fitness=Best_fitness;
save('best_SPBO_ROBL_T2pre_T2post.mat','SPBO_ROBL_T2pre_T2post_best_fitness');