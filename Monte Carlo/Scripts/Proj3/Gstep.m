function [thetaStep,lambdaStep] = Gstep(lambdaPrev, ni, tstep, Phi,d)
    thetaStep = fthetapost(lambdaPrev, Phi)
    lambdaStep = flambdapost(d, thetaStep, tstep, ni)
    
end

