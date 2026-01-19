%%%%%%%%%%% Reactive Power Dispatch %%%%%%%%%%%
%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%
%Iteration=200, Swarm size=50
clear
clc

iter=0;
iteration=200;
particlenumber=50;

%Inertia Weight
w_max=0.9;
w_min=0.4;
w_temp=w_max;
w_step=(w_max-w_min)/iteration;

%Maximum and minimum limit
vol_min=0.95;
vol_max=1.10;

%Load IEEE 14 bus data
[baseMVA, bus, gen, branch]=loadcase(case14);

%Initialization of Swarm & velocity
%Control variables: vg1 (1.06), vg2 (1.045), vg3 (1.01), vg6 (1.07), vg8(1.09)
%Control variables: tp1(4-7 0.978), tp2(4-9 0.969), tp3(5-6 0.932)
%Control variables: shunt9(19)

%Random 50*9 matrix
Swarm=[unifrnd(0.95,1.10,particlenumber,5), ...
unifrnd(0.975,1.025,particlenumber,3),unifrnd(0,20,particlenumber,1)];

%Initial velocity is set to 0
Velocity=zeros(particlenumber,9);

for i=1:particlenumber
     v1=Swarm(i,1); %v1
     bus(1,8)=v1; %Vm, 8 is voltage magnitude (p.u.)
     gen(1,6)=v1; %Vg, 6 is voltage magnitude setpoint (p.u.)
     v2=Swarm(i,2); %v2
     bus(2,8)=v2;
     gen(2,6)=v2;
    
     v3=Swarm(i,3); %v3
     bus(3,8)=v3;
     gen(3,6)=v3;
     v6=Swarm(i,4); %v6
     bus(6,8)=v6;
     gen(4,6)=v6;
     v8=Swarm(i,5); %v8
     bus(8,8)=v8;
     gen(5,6)=v8;
     branch(8,9)=Swarm(i,6); %tp1 4-7, 9 is tap position
     branch(9,9)=Swarm(i,7); %tp2 4-9
     branch(10,9)=Swarm(i,8); %tp3 5-6
     bus(9,6)=Swarm(i,9); %Shunt capacitor 9, 6 is BS
    
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

%%%%%%%%%%%%% Initialize Pbest and Gbest %%%%%%%%%%%%%

% update local best
for j=1:particlenumber
 Pbest(j,:)=Swarm(j,:);
 Val_Pbest(j)=Obj_fun_initial(j);
end

% update global best
[Val_Gbest,m]=min(Val_Pbest);
Gbest=Swarm(m,:);
Gbest_calc=repmat(Swarm(m,:),particlenumber,1);

%%%%%%%%%%%%%%%% PSO %%%%%%%%%%%%%%%%%%
losses_temp=zeros(1,particlenumber);
figure('NumberTitle', 'off', 'Name', 'Minimum Real Power Loss');
title('Minimum Real Power Loss');
ylabel('Real Power Loss (MW)');
xlabel('Iteration');
grid on;
hold on

for iter=1:iteration
     R1=rand(particlenumber,9);
     R2=rand(particlenumber,9);
     %2.05+2.05=4.1;
     %2/abs(2-4.1-sqrt(4.1*4.1-4*4.1))=0.729
     
     %velocity update
     Velocity=0.729*(w_temp*Velocity+2.05*R1.*(Pbest-Swarm)+2.05*R2.*(Gbest_calc-Swarm));
    
     % Set maximum velocity
     for v_iter=1:9
         if v_iter==9
             Outstep=Velocity(:,v_iter)>0.1;
             Velocity(find(Outstep),v_iter)=0.1;
             Outstep=Velocity(:,v_iter)<-0.1;
             Velocity(find(Outstep),v_iter)=-0.1;
         else
             Outstep=Velocity(:,v_iter)>0.003;
             Velocity(find(Outstep),v_iter)=0.003;
             Outstep=Velocity(:,v_iter)<-0.003;
             Velocity(find(Outstep),v_iter)=-0.003;
         end
     end
     
     % position update
     Swarm=Swarm+Velocity;
    
     for k=1:particlenumber
         v1=Swarm(k,1); %v1
         bus(1,8)=v1; %Vm, 8 is voltage magnitude (p.u.)
         gen(1,6)=v1; %Vg, 6 is voltage magnitude setpoint (p.u.)
         v2=Swarm(k,2); %v2
         bus(2,8)=v2;
         gen(2,6)=v2;
         v3=Swarm(k,3); %v3
         bus(3,8)=v3;
         gen(3,6)=v3;
         v6=Swarm(k,4); %v6
         bus(6,8)=v6;
         gen(4,6)=v6;
         v8=Swarm(k,5); %v8
         bus(8,8)=v8;
         gen(5,6)=v8;
         branch(8,9)=Swarm(k,6); %tp1 4-7, 9 is tap position
         branch(9,9)=Swarm(k,7); %tp
        
         branch(10,9)=Swarm(k,8); %tp3 5-6
         bus(9,6)=Swarm(k,9); %Shunt capacitor 10, 6 is BS

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
         losses_temp(k)=eval(['losses',num2str(k)]);
        
         Obj_fun_temp(k)=losses_temp(k)+penalty_bus_violation+penalty_gen_violation+penalty_brch_violation;
        
         if Obj_fun_temp(k)<Val_Pbest(k)
             losses(k)=losses_temp(k);
             Val_Pbest(k)=Obj_fun_temp(k);
             Pbest(k,:)=Swarm(k,:);
         end
     end
    
     [Val_Gbest_temp,n]=min(Val_Pbest);
     if Val_Gbest_temp<Val_Gbest
         Val_Gbest=Val_Gbest_temp;
         Gbest=Swarm(n,:);
         Gbest_calc=repmat(Swarm(n,:),particlenumber,1);
     end
    
     w_temp=w_temp-w_step;
     Val_Gbest_rec(iter)=Val_Gbest;
     plot(Val_Gbest_rec);
     drawnow;
end