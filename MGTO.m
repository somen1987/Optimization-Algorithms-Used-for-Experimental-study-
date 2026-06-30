%tic; 
       
%clear all;
%clc;

%fhd=str2func('cec17_func');

%load('extremum2017');
function out=MGTO(problem,params)
%function [Silverback_Score,Silverback,convergence_curve]=GTO(pop_size,max_iter,lower_bound,upper_bound,variables_no,fobj)
global fixedimg transformedimg;
CostFunction =problem.CostFunction;
pop_size=50;
lower_bound=[problem.TxMin problem.TyMin problem.RotMin];
upper_bound=[problem.TxMax problem.TyMax problem.RotMax];
variables_no=problem.nVar;

%target_err=1e-8;%when target error is greater thn 10^-8 then prog will terminate.
%TEST_RUN=51;
%fes=zeros(TEST_RUN,1);%initialize fes by zeros 
%best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=5000;
max_iter=params.MaxIt;


%Start_fun=input('Starting function#');
%End_fun=input('Ending function#');

%for func_num=Start_fun:1:End_fun

tic;
%for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
%filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',variables_no,func_num,test,variables_no);    

%load(filename);




% initialize Silverback
Silverback=[];
Silverback_Score=inf;

%Initialize the first random population of Gorilla
X=initialization(pop_size,variables_no,upper_bound,lower_bound);

fitcount=0;
convergence_curve=zeros(max_iter,1);

for i=1:pop_size   
    Pop_Fit(i,1)=feval(CostFunction,X(i,:)');
    fitcount=fitcount+1;
    if Pop_Fit(i)<Silverback_Score 
            Silverback_Score=Pop_Fit(i,1); 
            Silverback=X(i,:);
    end
end

cm(1)=0.7;
k=1;
for i=1:pop_size  
    
    
    for j=1:1:variables_no
        
     cm(k+1)=4.*cm(k).*(1-cm(k));
     k=k+1;
     CX(i,j)=lower_bound(j)+(upper_bound(j)-lower_bound(j)).*cm(k);
    end
    
    CPop_Fit(i,1)=feval(CostFunction,CX(i,:)');
    fitcount=fitcount+1;
    if CPop_Fit(i,1)<Silverback_Score 
            Silverback_Score=CPop_Fit(i,1); 
            Silverback=CX(i,:);
    end
end

   AllFitness=[Pop_Fit CPop_Fit];
   All_X=[X; CX];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:pop_size),:);
   Pop_Fit=sorted_fitness(1:pop_size);


for i=1:pop_size  
    OX(i,:)=upper_bound+lower_bound-X(i,:);
    OPop_Fit(i,1)=feval(CostFunction,OX(i,:)');
    fitcount=fitcount+1;
    if OPop_Fit(i,1)<Silverback_Score 
            Silverback_Score=OPop_Fit(i); 
            Silverback=OX(i,:);
    end
end
  
%   size(Pop_Fit)
%    size(OPop_Fit)
%    

   AllFitness=[Pop_Fit OPop_Fit'];
   All_X=[X; OX];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:pop_size),:);
   Pop_Fit=sorted_fitness(1:pop_size);
   


GX=X(:,:);
lb=ones(1,variables_no).*lower_bound; 
ub=ones(1,variables_no).*upper_bound; 



%err=abs(extremum(func_num) - Silverback_Score);

%%  Controlling parameter

p=0.03;
Beta=3;
w=0.8;
Pgj(1)=0.3;
%%Main loop
%for It=1:max_iter 
%err=Inf;
eta=1;
It=1;
while fitcount < FES  
    
    if  mod(It,50)==0
        beta=0.01+(0.9-0.01).*rand; 
        Pgj(It)=beta.* Pgj(It-1);
        eta=beta.*eta;
       %Pgj=sine_map(Pgj);
       %Pgj=tent_map(Pgj);
        %Pgj=cubic_map(Pgj);
        %Pgj=gaussian_map(Pgj);
        %Pgj=logistic_map(cm);
        %Pgj=0.4*normrnd(mu,sigma);
    elseif It>1
        
        Pgj(It)=Pgj(It-1);
        
    end    
    
    
    
    
    if rand < Pgj
    
for i=1:pop_size  
    OX(i,:)=upper_bound+lower_bound-X(i,:);
    OPop_Fit(i,1)=feval(CostFunction,OX(i,:)');
    fitcount=fitcount+1;
%     if OPop_Fit(i)<Silverback_Score 
%             Silverback_Score=OPop_Fit(i); 
%             Silverback=OX(i,:);
%     end
end


   AllFitness=[Pop_Fit OPop_Fit'];
   All_X=[X; OX];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:pop_size),:);
   Pop_Fit=sorted_fitness(1:pop_size);
   Silverback_Score=Pop_Fit(1);
   Silverback=X(1,:);
        
        
    else
    a=(cos(2*rand)+1)*(1-It/max_iter);
    C=a*(2*rand-1); 

%% Exploration:

    for i=1:pop_size
        if rand<p    
            GX(i,:) =(ub-lb)*rand+lb;
        else  
            if rand>=0.5
                Z = unifrnd(-a,a,1,variables_no);
                H=Z.*X(i,:);   
                GX(i,:)=(rand-a)*X(randi([1,pop_size]),:)+C.*H; 
            else   
                GX(i,:)=X(i,:)-C.*(C*(X(i,:)- GX(randi([1,pop_size]),:))+rand*(X(i,:)-GX(randi([1,pop_size]),:))); %ok ok 

            end
        end
    end       
       
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation 
    for i=1:pop_size
         New_Fit= feval(CostFunction,GX(i,:)');    
         fitcount=fitcount+1;
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
    
%% Exploitation:  
    for i=1:pop_size
       if a>=w  
            g=2^C;
            delta= (abs(mean(GX)).^g).^(1/g);
            GX(i,:)=C*delta.*(X(i,:)-Silverback)+X(i,:); 
       else
           
           if rand>=0.5
              h=randn(1,variables_no);
           else
              h=randn(1,1);
           end
           r1=rand; 
           GX(i,:)= Silverback-(Silverback*(2*r1-1)-X(i,:)*(2*r1-1)).*(Beta*h); 
           
       end
    end
   
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation    
    for i=1:pop_size
         New_Fit= feval(CostFunction,GX(i,:)');
         fitcount=fitcount+1;
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
    end
It=It+1;             
convergence_curve(It)=Silverback_Score;
%err=abs(extremum(func_num) - Silverback_Score);
%fprintf("In Iteration %d, best estimation of the global optimum is %4.4f \n ", It,Silverback_Score );
         
end 


fprintf('global best=%f\n',Silverback_Score);
 
% best_it(test,1)=err;
 
 
% fprintf('RUN %d is completed=================\n',test);
 
 %fes(test,1)=fitcount;
%end

%fprintf('Function#%d is completed=====>\n',func_num);

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
%filename='GTOOBL.xlsx';
%sheet=sprintf('%dD',variables_no);
%if func_num>=26
    %column=sprintf('A%s','A'+func_num-26);
%else
  %   column=sprintf('%s','A'+func_num);
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
out.BestSol=Silverback;
out.BestCost=Silverback_Score;
image_transformation(out.BestSol);

%imtool(transformedimg)
out.BestPcc=pcc(fixedimg,transformedimg);
end%func