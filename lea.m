function [VND, IC50, DMAX, VS] = lea(VTbase,VTblock,plasmaDrug,S,startVals)
% Estimates personal VND and global IC50 and max occupancy values for a
% dataset of baseline and block scans for N subjects using Likelihood 
% Estimation of Affinity (LEA). 
% 
% Inputs: 
%   VTbase:     kxN array of baseline VTs in k ROIs for N subjects 
%   VTblock:    kxN array of post-drug VTs in k ROIs for N subjects
%   plasmaDrug: 1xN array of measured plasma drug concentrations for N
%                   subjects
%   S:          kxk covariance matrix
%   startVals:  (N+2)x1 array of initial guesses for parameter values. 
%                   startVals(1:N) for VND
%                   startVals(N+1) for D_MAX (max attainable occupancy)
%                   startVals(N+2) for IC_50
%               Can be omitted, in which case the initial guesses will be
%               based on the Lassen Plot and dose-occupancy plot.
% Outputs: 
%   VND:        1xN array of estimated VND values
%   IC50:       estimated population IC_50
%   DMAX:       estimated population D_max 
%
%--------------------------------------------------------------------------
%                                              	Gjertrud Louise Laurell
%                                               Neurobiology Research Unit 
%                                               May 2022

if nargin < 5
    startVals = lea_getStartVals(VTbase,VTblock,plasmaDrug); 
end

N = size(VTbase,2);     % Number of subjects

% Constrain parameters, so that 0 < dmax < 1, 0 < vnd, and 0 < ic50
x0([1:N N+2]) = log(startVals([1:N N+2]));
x0(N+1) = log(startVals(N+1)/(1 - startVals(N+1))); 

options = optimset('MaxFunEvals',100000,'MaxIter',100000);

% Find maximum of log-likelihood function (l) by minimizing -l
x = fminsearch(@ll_LEA,x0,options,VTbase,VTblock,plasmaDrug,S); 

% Transform back parameter values
VND = exp(x(1:N)); 
DMAX = exp(x(N+1))/(1 + exp(x(N+1)));
IC50 = exp(x(N+2)); 

% Calculate VS 
O = 1 - DMAX*(plasmaDrug./(plasmaDrug + IC50)); 
VS = (VTbase - VND + O.*(VTblock - VND))./(1 + O.^2); 

function l = ll_LEA(x,vt1,vt2,cp,sigma)
% The LEA log-likelihood function. vt1 are the baseline VT values, vt2 are
%  the post-drug VT values, cp are the plasma drug concentrations, and 
% sigma is the covariance matrix

% Transform back parameter values 
vnd = exp(x(1:N));
dmax = exp(x(N+1))/(1+exp(x(N+1)));
ic50 = exp(x(N+2));

lj = zeros(N,1);    % Array of log-likelihood functions for each subject

for id = 1:N 
    oj = 1 - dmax*(cp(id)/(cp(id) + ic50)); 
    vsj = (vt1(:,id) - vnd(id) + oj*(vt2(:,id) - vnd(id)))/(1 + oj^2); 
    
    a = vt1(:,id) - vnd(id) - vsj; 
    base = a'*(sigma\a); 
    
    b = vt2(:,id) - vnd(id) - oj*vsj; 
    block = b'*(sigma\b);
    
    lj(id) = base + block; 
end

% Total log-likelihood is the sum of contributions from each subject
l = sum(lj); 

end

end