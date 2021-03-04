function [thetaStep,lambdaStep] = Gstep(lambdaPrev, thetaPrev, tau, tstep, Phi)
    thetaStep = fthetaPost(lambdaPrev, Phi);
    lambdaStep = flambdaPost(thetaStep, tstep, tau);
    
end

