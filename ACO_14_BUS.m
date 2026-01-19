%% Continuous Ant Colony Optimization (ACOR) for Rosenbrock Function
clc; clear; close all;

% Number of Decision Variables (x and y)
nVar = 9;               

VarMin = [0.95 0.95 0.95 0.95 0.95 0.975 0.975 0.975 0.0]; 
VarMax = [1.10 1.10 1.10 1.10 1.10 1.025 1.025 1.025 20];

% ACO Parameters
MaxIt = 100;            % Maximum Number of Iterations
nAnt = 50;              % Number of Ants
nArchive = 60;          % Size of the Solution Archive
q = 0.1;                % Intensification Factor (Selection pressure)
zeta = 1.1;               % Deviation-distance Ratio (Exploration)

% Initialization
% Create empty individual structure
template.Position = [];
template.Cost = [];

% Create Archive
archive = repmat(template, nArchive, 1);

% Initialize Archive with random solutions
for i = 1:nArchive
    archive(i).Position = unifrnd(VarMin, VarMax);
    %%%% < function Call > %%%%%
    archive(i).Cost = LossFunction(archive(i).Position);
end

% Sort Archive
[~, SortOrder] = sort([archive.Cost]);
archive = archive(SortOrder);

% Best Solution ever found
BestSol = archive(1);

% Array to hold Best Cost values
BestCost = zeros(MaxIt, 1);


% --- Setup Dynamic Plotting ---
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]);

% Subplot 1: Convergence
subplot(1, 2, 1);
hConv = plot(NaN, NaN, 'r', 'LineWidth', 2);
xlabel('Iteration'); ylabel('Best Cost');
title('Convergence Curve');
grid on;

% Subplot 2: Population Distribution (Dimensions 1 and 2)
subplot(1, 2, 2);
hPop = plot(NaN, NaN, 'b.', 'MarkerSize', 10);
hold on;
hBest = plot(NaN, NaN, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xlabel('x_1'); ylabel('x_2');
title('Ant Population (Dimensions 1 & 2)');
axis([-2 2 -2 2]);
grid on;
% ------------------------------


% 4. Main Loop
for it = 1:MaxIt

    % Calculate Weights (Pheromone intensity)
    w = 1/(q*nArchive*sqrt(2*pi)) * exp(-0.5*(((1:nArchive)-1)/(q*nArchive)).^2);
    p = w/sum(w); % Probabilities for choosing a solution from archive

    % New Population
    pop = repmat(template, nAnt, 1);

    for i = 1:nAnt
        % Select a guide (Roulette wheel selection)
        l = find(rand <= cumsum(p), 1, 'first');

        % Generate New Position
        pop(i).Position = zeros(1, nVar);
        for j = 1:nVar
            % Calculate standard deviation
            sigma = 0;
            for r = 1:nArchive
                sigma = sigma + abs(archive(r).Position(j) - archive(l).Position(j));
            end
            sigma = zeta * sigma / (nArchive - 1);

            % Sample from Gaussian Distribution
            pop(i).Position(j) = archive(l).Position(j) + sigma * randn;

            % Boundary Constraints
            pop(i).Position(j) = max(pop(i).Position(j), VarMin(j));
            pop(i).Position(j) = min(pop(i).Position(j), VarMax(j));
        end

        % Evaluation
        %%%%%% < function Call > %%%%%%
        pop(i).Cost = LossFunction(pop(i).Position);
    end

    % Update Archive
    archive = [archive; pop];
    [~, SortOrder] = sort([archive.Cost]);
    archive = archive(SortOrder);
    archive = archive(1:nArchive); % Keep only the best

    % Update Best Solution
    BestSol = archive(1);
    BestCost(it) = BestSol.Cost;

    % % Display iteration info
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);


    % --- Update Dynamic Plots ---
    set(hConv, 'XData', 1:it, 'YData', BestCost(1:it));

    % Get all ant positions for the second plot
    allPos = vertcat(pop.Position);
    set(hPop, 'XData', allPos(:,1), 'YData', allPos(:,2));
    set(hBest, 'XData', BestSol.Position(1), 'YData', BestSol.Position(2));

    drawnow; % Forces MATLAB to update the figure window
    % ----------------------------

end

%% 5. Results
fprintf('\n--- Optimization Results ---\n');
disp('Best Position Found:');
disp(BestSol.Position);
fprintf('Minimum Value: %0.6e\n', BestSol.Cost);




function [Obj_fun] = LossFunction(x)
    [baseMVA, bus, gen, branch]=loadcase(case14);
    %Maximum and minimum limit
    vol_min=0.95;
    vol_max=1.10;

         bus(1,8)=x(1); %Vm, 8 is voltage magnitude (p.u.)
         gen(1,6)=x(1); %Vg, 6 is voltage magnitude setpoint (p.u.)
         
         bus(2,8)=x(2);
         gen(2,6)=x(2);

         bus(3,8)=x(3); 
         gen(3,6)=x(3); 
         
         bus(6,8)=x(4); 
         gen(4,6)=x(4); 
         
         bus(8,8)=x(5); 
         gen(5,6)=x(5); 
         
         branch(8,9)=x(6); %tp1 4-7, 9 is tap position
         
         branch(9,9)=x(7); %tp2 4-9
         
         branch(10,9)=x(8); %tp3 5-6
         
         bus(9,6)=x(9); %Shunt capacitor 9, 6 is BS

         eval(['savecase (''case14_test.mat'', baseMVA, bus, gen, branch)']);
         eval(['results','=runpf(''case14_test.mat'')']);
         eval(['losses','=sum(real(get_losses(results)))']);
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
         losses=eval(['losses']);

         Obj_fun=losses+penalty_bus_violation+penalty_gen_violation+penalty_brch_violation;
end
