clear all 
close all 
clc

%% ME6103 Hw 2 KKT AND STOCHASTIC OPT TECHNIQUES SHUBH RAVAL

%% Problem 1 KKT and simulated anealing 

%% A) 
    % Assume Water Storage Tank to be a Cylinder 
    % Minimize the overall surface area of the tank: SA =  2*pi*r +
        % 2*pi*r*h [objective function]

    % Subject to: 
        % Volume of tank (V = pi * r^2 * h) >= 200m^3 [inequality
        % constraint]
        % ( height (h) / 2* radius (r) ) = 2 [equality constraint]

%% B) 
    % Formulate the Lagrangian: 
    % Thus the objective function can be written as: 
        % SA = 2*pi*r*h + 2*pi*r^2

    % The equality constraint is given by the ratio of height to
    % diameter (here written as 2*r) which is defined as: 
    % (h-4r) = 0 -> lamda*(h-4*r)

    % The inequality constraint finally is given by the Volume
    % needing to be greater than 200 m^3 which is defined as: 
    % V = pi * r^2 * h >= 200 -> mu * (200 - pi * r^2 * h)
    
    % This gives the full Lagrangian as:
    % L = (2*pi*r*h + 2*pi*r^2) + lamda*(h - 4*r) + mu*(200 - (pi*r^2*h))
    
%% C 
    % Derive the KKT conditions for optimality: 

    % KKT conditions to meet for optimiality is then defined as: 
    % dL/dx1 = df/dx1 + mu1*dg/dx1 + lamda1*dh/dx1 = 0
    % dL/dx2 = df/dx2 + mu1*dg/dx2 + lamda1*dh/dx2 = 0
    
    % For this problem: 
    % dL/dr = 2*pi*h + 4*pi*r + lamda*(-4) + mu*(-2*pi*r*h)
    % dl/dh = 2*pi*r + lamda + mu*(-pi*r^2)

%% D
    % Determine if KKT conditions are met for the following points: 

    %% i) r = 2.515 m , h = 10.06 m 

        % dL/dr = 0 = 2*pi*10.06 + 4*pi*2.515 + l*(-4) + mu*(-2*pi*2.515*10.06)
        % dL/dh = 0 = 2*pi*2.515 + l + mu*(-pi*(2.515)^2)

        % solving the system gives lamda = -2.6337 and mu = .662691

        % non negativity passes as mu >= 0

        % g = 0.0949 fails feasibility 

        % Complementarity fails as mu*g = 0.0629

        % h = 10.06 - 2.515*4 = 0 passes feasibility

        % So this very closely fails the KKT conditions but considering
        % that the margin for failure was 0.0629 since the volume was
        % 0.0949 m^3 smaller than desired means that practically this could
        % be a potential option for exploring slightly more to reach
        % optimality 


    %% ii) r = 1.25 m , h = 5.0 m 

        % dL/dr = 2*pi*5.0 + 4*pi*1.25 -4*l + mu*(-2*pi*1.25*5.0)
        % dL/dh = 0 = 2*pi*1.25 + l + mu*(-pi*(1.25)^2)

        % solving the system gives lamda = -1.309 and mu = 1.3333

        % non negativity passes as mu >= 0

        % g = 175.45 fails feasibility 

        % Complementarity fails as mu*g = 175.4563

        % h = 5 - 1.25*4 = 0 passes feasibility 

        % This dramatically fails the volume feasibility and
        % complementarity since the volume found was about 25 m^3 compared
        % to the 200 m^3 required


    %% iii) r = 4 m , h = 16 m 

        % dL/dr = 2*pi*16 + 4*pi*4 -4*l + mu*(-2*pi*4*16)
        % dL/dh = 0 = 2*pi*4 + l + mu*(-pi*(4)^2)

        % solving the system gives lamda = -4.1888 and mu = .41667

        % non negativity passes as mu >= 0

        % g = -604.24 passes feasibility 

        % mu*g = -251.7718 fails Complementarity 

        % h = 16 - 4*4 = 0 so passes feasibility

        % This also fails the KKT conditions even though g pass feasibility with
        % a volume of about 800 m^3 compared to the needed 200^3 this
        % however causes Complementarity to fail since in no way can mu
        % bring this value close to or at 0
  
 %% 2 Solve the test problem via a manual genetic algorithm 

    f =  @(x) sin(x) + (0.05*x.^2) + 1; 

    %% Parameters: 

    popSize = 10; 
    numGen = 5; 
    crossProb = 0.6;
    mutProb = 0.05; 
    numBits = 6;
    xMin = -7;
    xMax = 7;

    %% Starting Population Initialized as Binary 
    binaryPop = randi([0,1], popSize, numBits); 
    s = rng; 
    
    %% Convert to Decimal and Scale to the xMin and xMax range
    convert_scale_binary = @(bin) xMin + (bin2dec(num2str(bin)) / (2^numBits - 1)) * (xMax - xMin); 

    %% Best Sol
    bestSolutions = zeros(numGen,2);

    figure; 
    hold on; 
    x = linspace(xMin, xMax, 100); 
    plot(x, f(x), 'b-');


    gen = 1; 
    legendLabels = {'Function'};  
    while gen <= numGen
          decodedPop = arrayfun(@(row) convert_scale_binary(binaryPop(row, :)), 1:popSize); 
          fitness = f(decodedPop); 

          [~,sortedIndices] = sort(fitness); 
          binaryPop = binaryPop(sortedIndices, :);
          decodedPop = decodedPop(sortedIndices);
          fitness = fitness(sortedIndices);

          %bestSol(gen, :) = [decodedPop(1), fitness(1)];
          scatter(decodedPop, fitness, 'filled');

          title('Genetic Algorithm');
          pause(1); 

          fprintf('\nGeneration %d\n', gen);
          fprintf('| %10s | %20s | %10s | %10s |\n', 'Fitness', 'Binary', 'Decimal', 'Generation');
          fprintf('|------------|----------------------|------------|------------|\n');
        
          for i = 1:popSize
            fprintf('| %10.4f | %20s | %10.4f | %10d |\n', fitness(i), num2str(binaryPop(i, :)), decodedPop(i), gen);
          end
          legendLabels{length(legendLabels) + 1} = sprintf('Generation %d', gen);
          %legend show; 
          % Elitism
          newPop = binaryPop(1:2,:);

          %Roulette
          totalFit = sum(1 ./ fitness); %inverse since minimize
          selectionProbs = (1 ./ fitness) / totalFit; 
          cumProbs = cumsum(selectionProbs); 

          %Make the kids
          while size(newPop, 1) < popSize
              r1 = rand();
              r2 = rand();
              parent1 = binaryPop(find(cumProbs >= r1, 1), :);
              parent2 = binaryPop(find(cumProbs >= r2, 1), :);

              %Cross
              if rand() < crossProb
                  crossoverPoint = randi(numBits); 
                  kid1 = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
                  kid2 = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
              else
                  kid1 = parent1;
                  kid2 = parent2;
              end

              % Mutation

              for i = 1:numBits
                  if rand() < mutProb, kid1(i) = ~kid1(i); end
                  if rand() < mutProb, kid2(i) = ~kid2(i); end
              end

              % new pop

              newPop = [newPop; kid1; kid2];
              if size(newPop,1)>popSize
                  newPop = newPop(1:popSize, :);
              end
          end
              binaryPop = newPop; 
              gen = gen+1; 
    end 


    legend(legendLabels,'Location', 'best');
    hold off;
              






