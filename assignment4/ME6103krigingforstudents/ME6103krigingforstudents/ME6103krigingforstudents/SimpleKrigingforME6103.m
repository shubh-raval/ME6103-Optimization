% Simple kriging fit to a sine function
% 
%===============

x1 = [0.6283; 1.2566; 1.8850; 2.5133; 3.1416; 3.7699; 4.3982; 5.0265; 5.6549; 6.2832];
x2 = [ 3.1416; 1.5708; 4.7124; 0.7854; 3.9270; 2.3562; 5.4978; 0.3927; 3.5343; 1.9635];
xTrain = [x1 x2]; 
% xTrain is the set of design variable values for training the kriging model.  
%Each row of the xTrain matrix should contain the x1 and x2 values for a unique input point.  
%The number of rows should correspond to the number of training points.
yTrain = sin(x1) + sin(x2); 
%yTrain is the set of responses for training the kriging model. 
%It should be an n x 1 vector, where n = number of training points.

% set the kriging model parameters (See the dace.pdf file for info, Sec 5
% and 6 for more info)
regr = @regpoly0;
corr = @corrgauss;
theta0 = 0.1;
lob = 0.000001;
upb = 10;

% build the dacefit model
[dmodel, perf] = dacefit(xTrain, yTrain, regr, corr, theta0, lob, upb);

% test some new points
xTest =  [1.2566 3.1416; 2.5133 1.5708; 3.7699 4.7124; 5.0265 0.7854; 6.2832 3.9270];
yTest = predictor(xTest, dmodel);