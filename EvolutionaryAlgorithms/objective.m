%% Define the objective function 
function obj = objective(p,act,refTheta,time)
    % Model predicted output
    thetaP = simulation(p,act,refTheta,time);  
    % Calculate the Objective
    % mean square error
    % output the objective value
    error = refTheta - thetaP;
    obj = sum(error.^2)/length(thetaP);     
 end