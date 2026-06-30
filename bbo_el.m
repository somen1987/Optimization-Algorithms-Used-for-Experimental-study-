%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA113
% Project Title: Biogeography-Based Optimization (BBO) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%


%% Problem Definition
function out=bbo_el(problem,params)
global fixedimg transformedimg
CostFunction =problem.CostFunction;        % Cost Function

nVar=3;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=[problem.TxMin problem.TyMin problem.RotMin];         % Decision Variables Lower Bound
VarMax=[problem.TxMax problem.TyMax problem.RotMax];         % Decision Variables Upper Bound

%% BBO Parameters

MaxIt=100;          % Maximum Number of Iterations

nPop=50;            % Number of Habitats (Population Size)

KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

nNew=nPop-nKeep;                % Number of New Habitats

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates
lambda=1-mu;                    % Immigration Rates

alpha=0.7;

pMutation=0.1;

sigma=0.02*(VarMax-VarMin);

%% Initialization

% Empty Habitat
habitat.Position=[];
habitat.Cost=[];

% Create Habitats Array
pop=repmat(habitat,nPop,1);
VMax=0.3.*(VarMax-VarMin);
VMin=VMax;
c1=2.05;
c2=c1;
% Initialize Habitats
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Velocity=unifrnd(VMin,VMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    pop(i).pbest=pop(i).Cost;
end

% Sort Population
[~, SortOrder]=sort([pop.Cost]);
pop=pop(SortOrder);

% Best Solution Ever Found
BestSol=pop(1);
gbest=BestSol.Position;
gb=BestSol.Cost;
% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop

for it=1:MaxIt
    
    newpop=pop;
    for i=1:nPop
        for k=1:nVar
            % Migration
            if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;
                EP=EP/sum(EP);
                
                % Select Source Habitat
                j=RouletteWheelSelection(EP);
                
                % Migration
                
%                 newpop(i).Position(k)=pop(i).Position(k) ...
%                     +alpha*(pop(j).Position(k)-pop(i).Position(k));
                
                r1=round(nPop*rand);
                
                while (i~=r1) r1=round(nPop*rand); end
                
                r2=round(nPop*rand);
                    
                while (i~=r2 && r1~=r2) r2=round(nPop*rand); end
                
                r3=round(nPop*rand);
                    
                while (i~=r3 && r3~=r1 && r3~=r2) r3=round(nPop*rand); end
                
                newpop(i).Position(k)=pop(r1).Position(k) ...
                     +alpha*(pop(r2).Position(k)-pop(r3).Position(k));
                
            end
            
            % Mutation
            if rand<=pMutation
                newpop(i).Position=newpop(i).Position(k)+sigma*randn;
            end
        end
        
        if (sum(abs(newpop(i).Position-pop(i).Position))==0)
            %disp('no-update');
            pop(i).Velocity=newpop(i).Velocity+c1.*rand(1,nVar).*(pop(i).pbest-pop(i).Position)...
                +c2.*rand(1,nVar).*(gbest-pop(i).Position);
            pop(i).Velocity=max(min(pop(i).Velocity,VMax),VMin);
            newpop(i).Position=pop(i).Position +  pop(i).Velocity; 
        end    
        
        % Apply Lower and Upper Bound Limits
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation
        newpop(i).Cost=CostFunction(newpop(i).Position);
        
        
        if  newpop(i).Cost < pop(i).pbest
            newpop(i).pbest=newpop(i).Position;
            pop(i).pbest=newpop(i).Cost;
        end
        
        if  newpop(i).Cost < gb
            gbest=newpop(i).Position;
            gb=newpop(i).Cost;
        end
        
        ConvergenceCurve(i)=newpop(i).Cost;
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population
    pop=[pop(1:nKeep)
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    gbest=BestSol.Position;
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
out.BestSol=gbest;
out.BestCost=BestCost(it);
image_transformation(out.BestSol);
out.ConvergenceCurve(it)=BestCost(it);
%imtool(transformedimg)
out.BestPcc=pcc(fixedimg,transformedimg);
 image_transformation(out.BestSol);

%imtool(transformedimg)
out.Bestssim=ssim(fixedimg,transformedimg); 
image_transformation(out.BestSol);

out.Bestncc=NCC(fixedimg,transformedimg);
image_transformation(out.BestSol);

out.Bestecc=ECC(fixedimg,transformedimg);
image_transformation(out.BestSol);

end







