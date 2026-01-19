%%%% FPA Based Reactive Power Dispatch Without New DGs%%%% %%%%%%%%%%%%%%%Initialization %%%%%%%%%%%%%% %Iteration=200, %Particle size=50 

clear;
clc;
 
iter=0;
iteration=200; 
particlenumber=50;
probability = 0.8; %
 
d=9;           % Dimension of the search variables
 
%Maximum and minimum limit 
vol_min=0.95;
vol_max=1.10;
Lb = [0.95 0.95 0.95 0.95 0.95 0.975 0.975 0.975 0];
Ub = [1.10 1.10 1.10 1.10 1.10 1.025 1.025 1.025 20]; 
 
%Load IEEE 14 bus data
[baseMVA, bus, gen, branch]=loadcase(case14);
 
%Initialization of pollen & velocity
%Control variables: vg1 (1.06), vg2 (1.045), vg3 (1.01), vg6 (1.07), vg8(1.09)
%Control variables: tp1(4-7 0.978), tp2(4-9 0.969), tp3(5-6 0.932)
%Control variables: shunt9(19)
 
%Random 50*9 matrix 
pollen=[unifrnd(0.95,1.10,particlenumber,5), ...
unifrnd(0.975,1.025,particlenumber,3),unifrnd(0,20,particlenumber,1)];
 
 
    for i=1:particlenumber 
        v1=pollen(i,1);  %v1
        bus(1,8)=v1;    %Vm, 8 is voltage magnitude (p.u.) 
        gen(1,6)=v1;    %Vg, 6 is voltage magnitude setpoint (p.u.) 
        v2=pollen(i,2);  %v2
        bus(2,8)=v2;
        gen(2,6)=v2;
 
        v3=pollen(i,3);  %v3 
        bus(3,8)=v3;
        gen(3,6)=v3;   
        v6=pollen(i,4);  %v6 
        bus(6,8)=v6;
        gen(4,6)=v6; 
        v8=pollen(i,5);  %v8 
        bus(8,8)=v8;
        gen(5,6)=v8;
 
        branch(8,9)=pollen(i,6); %tp1 4-7, 9 is tap position 
        branch(9,9)=pollen(i,7); %tp2 4-9
        branch(10,9)=pollen(i,8); %tp3 5-6 
        bus(9,6)=pollen(i,9);    %Shunt capacitor 9, 6 is BS
        eval(['savecase (''case14_test' num2str(i) '.mat'', baseMVA, bus, gen, branch)']); 
        eval(['results',num2str(i),'=runpf(''case14_test', num2str(i) '.mat'')']);
        eval(['losses',num2str(i),'=sum(real(get_losses(results',num2str(i),')))']);
 
        %Penalty for bus voltage violation 
        bus_inf=bus(:,8);
        for bus_num=1:14
            if bus_inf(bus_num)>vol_max 
                penalty_bus(bus_num)=10000*(bus_inf(bus_num)-vol_max)^2;
            elseif bus_inf(bus_num)<vol_min 
                penalty_bus(bus_num)=10000*(bus_inf(bus_num)-vol_min)^2;
            else
                penalty_bus(bus_num)=0; 
            end
        end
        penalty_bus_violation=sum(penalty_bus);
 
        %Penalty for reactive generation violation 
        gen_inf=gen(:,3);
        for gen_num=2:5
            if gen_inf(gen_num)>gen(gen_num,4) 
                penalty_gen(gen_num)=1000*(gen_inf(gen_num)-gen(gen_num,4))^2;
            elseif gen_inf(gen_num)<gen(gen_num,5) 
                penalty_gen(gen_num)=1000*(gen_inf(gen_num)-gen(gen_num,5))^2;
            else
                penalty_gen(gen_num)=0; 
            end
 
        end
        penalty_gen_violation=sum(penalty_gen);
 
        %Penalty for tap position violation 
        brch_inf=[branch(8,9); branch(9,9); branch(10,9)]; 
        for brch_num=1:3
            if brch_inf(brch_num)>1.025 
                penalty_brch(brch_num)=10000*(brch_inf(brch_num)-1.025)^2;
            elseif brch_inf(brch_num)<0.975 
                penalty_brch(brch_num)=10000*(brch_inf(brch_num)-0.975)^2;
            else
                penalty_brch(brch_num)=0; 
            end
        end
        penalty_brch_violation=sum(penalty_brch);
 
        %Penalty function 
        losses(i)=eval(['losses',num2str(i)]);
 
    Obj_fun_initial(i)=losses(i)+penalty_bus_violation+penalty_gen_violation+penalty_brch_violation; 
    
    end
    
    for j=1:particlenumber
        FP_best(j,:)= pollen(j,:); 
        Val_FP_best(j)=Obj_fun_initial(j);
    end
 
        [FG_best,m]=min(Val_FP_best); 
        FGbest=pollen(m,:); 
        FGbest_calc=repmat(pollen(m,:),particlenumber,1);
        losses_temp=zeros(1,particlenumber);
        
        
figure('NumberTitle', 'off', 'Name', 'Minimum Real Power Loss by FPA'); 
title('Minimum Real Power Loss');
ylabel('Real Power Loss (MW)'); 
xlabel('Iteration');
grid on; 
hold on
 
 
S = pollen;
 
 
% Start the main iterations -- Flower Pollination Algorithm 
for iter=1:iteration
    
        % Loop over all bats/solutions
        for k=1:particlenumber,
            v1=pollen(k,1); %v1
            bus(1,8)=v1;    %Vm, 8 is voltage magnitude (p.u.) 
            gen(1,6)=v1;    %Vg, 6 is voltage magnitude setpoint (p.u.) 
            v2=pollen(k,2);     %v2
            bus(2,8)=v2;
            gen(2,6)=v2; 
            v3=pollen(k,3); %v3 
            bus(3,8)=v3;
            gen(3,6)=v3; 
            v6=pollen(k,4); %v6 
            bus(6,8)=v6;
            gen(4,6)=v6; 
            v8=pollen(k,5); %v8 
            bus(8,8)=v8;
            gen(5,6)=v8;
 
            branch(8,9)=pollen(k,6); %tp1 4-7, 9 is tap position 
            branch(9,9)=pollen(k,7); %tp2 4-9
 
 
            branch(10,9)=pollen(k,8); %tp3 5-6 
            bus(9,6)=pollen(k,9);   %Shunt capacitor 10, 6 is BS
            
            eval(['savecase (''case14_test' num2str(k) '.mat'', baseMVA, bus, gen, branch)']); 
            eval(['results',num2str(k),'=runpf(''case14_test', num2str(k) '.mat'')']);
            eval(['losses',num2str(k),'=sum(real(get_losses(results',num2str(k),')))']);
 
            %Penalty for bus voltage violation 
            bus_inf=bus(:,8);
            for bus_num=1:14
                if bus_inf(bus_num)>vol_max 
                    penalty_bus(bus_num)=10000*(bus_inf(bus_num)-vol_max)^2;
                elseif bus_inf(bus_num)<vol_min 
                   penalty_bus(bus_num)=10000*(bus_inf(bus_num)-vol_min)^2;
                else
                    penalty_bus(bus_num)=0; 
                end
            end
            penalty_bus_violation=sum(penalty_bus);
 
            %Penalty for reactive generation violation gen_inf=gen(:,3);
            for gen_num=2:5
                if gen_inf(gen_num)>gen(gen_num,4) 
                    penalty_gen(gen_num)=1000*(gen_inf(gen_num)-gen(gen_num,4))^2;
                elseif gen_inf(gen_num)<gen(gen_num,5) 
                    penalty_gen(gen_num)=1000*(gen_inf(gen_num)-gen(gen_num,5))^2;
                else
                    penalty_gen(gen_num)=0; 
                end
            end
            penalty_gen_violation=sum(penalty_gen);
 
            %Penalty for tap position violation 
            brch_inf=[branch(8,9); branch(9,9); branch(10,9)]; 
            for brch_num=1:3
                if brch_inf(brch_num)>1.025 
                    penalty_brch(brch_num)=10000*(brch_inf(brch_num)-1.025)^2;
                elseif brch_inf(brch_num)<0.975 
                    penalty_brch(brch_num)=10000*(brch_inf(brch_num)-0.975)^2;
                else
                    penalty_brch(brch_num)=0;
                end
            end
            penalty_brch_violation=sum(penalty_brch);
 
            %Penalty function 
            losses_temp(k)=eval(['losses',num2str(k)]);
        
%          if k<2
%             old_obj_fun_temp = zeros(1,particlenumber);
%             old_obj_fun_temp = inf;
%          end
%          
%         Obj_fun_temp(k) =  old_obj_fun_temp;
%         old_obj_fun_temp = Obj_fun_temp(k);
        
             
        Obj_fun_temp(k)=losses_temp(k)+penalty_bus_violation+penalty_gen_violation+penalty_brch_violation;
            
            if Obj_fun_temp(k)<Val_FP_best(k) 
                %losses(k)=losses_temp(k); 
                Val_FP_best(k)=Obj_fun_temp(k); 
                FP_best(k,:)= pollen(k,:);
            end
   
            
          if rand>probability     % Probability (p) is checked after drawing a rand
                L=Levy(d);     % Draw random numbers from a Levy distribution
          % The main search mechanism via flower pollination
                dS=L.*(pollen(k,:)-FGbest);      % Caclulate the step increments
                S(k,:)=pollen(k,:)+dS;         % Update the new solutions
          % Check new solutions to satisfy the simple limits/bounds
                S(k,:)=simplebounds(S(k,:),Lb,Ub);
          % Another search move via local pollenation of neighbor flowers 
          else
              epsilon=rand;
              % Find random pollen/flowers in the neighbourhood
              JK=randperm(particlenumber);
              
              S(k,:)=S(k,:)+epsilon*(pollen(JK(1),:)-pollen(JK(2),:));
              % Check if the simple limits/bounds are met
              S(k,:)=simplebounds(S(k,:),Lb,Ub);
          end
          
          % Evaluate the objective values of the new solutions
           
           Fnew = Obj_fun_temp(k);
            if Fnew<= Obj_fun_temp(k),
                pollen(k,:)=S(k,:);
                Obj_fun_temp(k)=Fnew;
            end
           
             [FGbest_temp,n]=min(Val_FP_best); 
            if FGbest_temp<FG_best
                FG_best = FGbest_temp; 
                FGbest=pollen(n,:);
                FGbest_calc=repmat(pollen(n,:),particlenumber,1); 
            end
         end
   Val_FGbest_rec(iter)=FG_best; 
   plot(Val_FGbest_rec);
   drawnow;     
 end
        % Display results every 100 iterations
        if ~mod(iter,1)
        disp(strcat('Iteration iter=',num2str(iter)));
        % best
        FGbest
        end     
  % End of the main iterations
%% Output/display the final results
disp(['Total number of evaluations: ',num2str(iteration*particlenumber)]);
disp(['Best solution=',num2str(FGbest)]);
disp(['FG_best =',num2str(FG_best)]);
  
% Application of simple lower bounds and upper bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
end
% Draw n samples for the Levy flight from the Levy distribution
function L=Levy(d)
% For details of the Levy flights, see Chapter 11 of the following book:
% Xin-She Yang, Nature-Inspired Optimization Algorithms, Elsevier, (2014).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% Mantegna's algorithm for Levy random numbers
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;             % Final Levy steps
end
