clear all close all clc 

% % Test Binary Decoding and Fitness Calculation
% numBits = 6;      % Number of bits
% xMin = -7;        % Minimum x value
% xMax = 7;         % Maximum x value
% 
% % Choose a random binary string (just for testing)
% binaryStr = [0 1 1 0 1 0];  % Example binary string
% 
% % Convert binary string to decimal (manually)
% binStr = num2str(binaryStr);  % Convert binary vector to string
% binDec = bin2dec(binStr);     % Convert to decimal
% decimalValue = xMin + (binDec / (2^numBits - 1)) * (xMax - xMin);
% 
% % Define the fitness function
% f = @(x) sin(x) + (0.05 * x.^2) + 1;
% 
% % Calculate fitness manually
% manualFitness = f(decimalValue);
% 
% % Display results
% fprintf('Binary: %s\n', binStr);
% fprintf('Decoded Decimal: %.4f\n', decimalValue);
% fprintf('Fitness from Function: %.4f\n', manualFitness);


% Test Fitness Calculation
f = @(x) sin(x) + (0.05 * x.^2) + 1;

% Known values to test the fitness function
testValues = [-1.2222, -3, 0, 3, 7];  % Known x values within the range

% Calculate fitness for each of these test values
manualFitnessResults = arrayfun(f, testValues);

% Display the results
fprintf('Test Values and Corresponding Fitness:\n');
for i = 1:length(testValues)
    fprintf('x = %.4f, Fitness = %.4f\n', testValues(i), manualFitnessResults(i));
end