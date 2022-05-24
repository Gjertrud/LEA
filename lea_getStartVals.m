function x0 = lea_getStartVals(vtBaseline,vtBlock,plasmaDrug)
% Returns initail guesses or VND, max occupancy and IC50. 
% For VND, the Lassen Plot VND is used, as long as it positive, and smaller
% than the lowest baseline VT value. If not, half of the smallest baseline
% VT is used. 
%
% Inputs: 
%   vtBaseline: kxN array of baseline VTs in k ROIs for N subjects 
%   vtBlock:    kxN array of post-drug VTs in k ROIs for N subjects
%   plasmaDrug: 1xN array of measured plasma drug concentrations for N
%                   subjects
%   
% Output: 
%   x0: (N+2)x1 array of initial guesses for parameter values. 
%     	x0(1:N) for VND
%      	x0(N+1) for D_MAX (max attainable occupancy)
%      	x0(N+2) for IC_50
%
%--------------------------------------------------------------------------
%                                              	Gjertrud Louise Laurell
%                                               Neurobiology Research Unit 
%                                               May 2022

N = size(vtBaseline,2); 

% Calculate subject VND and occupancies with the Lassen plot
vndLassen = zeros(1,N);
vndGuess = zeros(1,N);
deltaLassen = zeros(1,N);
for n = 1:N 
    [vndLassen(n), deltaLassen(n)] = ...
        doLassenPlot(vtBaseline(:,n),vtBlock(:,n));
    
    % Initial guess for VND is based on the result of the Lassen plot
    if vndLassen(n) > 0 && vndLassen(n) < min(vtBaseline(:,n))
        vndGuess(n) = vndLassen(n); 
    else
        vndGuess(n) = min(vtBaseline(:,n))/2; 
    end
end

% Calculate Dmax and IC50 with dose-occupancy plot 
[dmaxLassen, ic50Lassen] = ...
    doseOccupancyResponseCurve(deltaLassen,plasmaDrug); 

% Set initail guess of Dmax based on result from dose-occupancy plot
if dmaxLassen > 0 && dmaxLassen < 1
    dmaxGuess = dmaxLassen; 
else
    dmaxGuess = 0.5; 
end

% Set initial guess of IC50 based on results from dose-occupancy plot
if ic50Lassen > 0 && ic50Lassen < max(plasmaDrug)
    ic50Guess = ic50Lassen; 
else
    ic50Guess = mean(plasmaDrug); 
end

x0 = [vndGuess'; dmaxGuess; ic50Guess]; 
end

% Functions 
function [VND,Delta] = doLassenPlot(VT_base,VT_block)
% Estimates VND and Occupancy with Lassen plot
%   VT_base     kx1 array of VT values at baseline
%   VT_block    kx1 array of VT values after block

Lassen = polyfit(VT_base,(VT_base-VT_block),1); 
Delta = Lassen(1); 
VND = -Lassen(2)/Delta; 
end

function [Dmax, IC50] = doseOccupancyResponseCurve(Delta,Cp)
% Calculates max occupancy (Dmax) and Kd given a set of corresponding
% occupansies (Delta) and plasma drug concentrations (Cp) by fitting a dose
% occupancy response curve
% Delta = Dmax*(Cp/(Cp + Kd))
if any(isnan(Delta))
    nanID = find(isnan(Delta));
    Cp(nanID) = []; Delta(nanID) = [];
end
if nnz(Cp) < 2
    Dmax = NaN;     IC50 = NaN;
else
    if nnz(Cp) ~= length(Cp)
        zID = find(Cp==0);
        Cp(zID) = []; Delta(zID) = [];
    end
    d0 = mean([max(Delta) 1]);
    k0 = mean(Cp);
    doseFun = @(b,x)((b(1)*x)./(b(2) + x));
    options = optimoptions('lsqcurvefit','Display','off');
    beta = lsqcurvefit(doseFun,[d0 k0],Cp,Delta,[0 0],[inf inf],options);

    Dmax = beta(1);     IC50 = beta(2);
end
end