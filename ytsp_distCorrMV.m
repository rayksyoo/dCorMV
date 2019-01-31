function [dCor, dCov, dVar, stats] = ytsp_distCorrMV(X, optCent)
% YTSP_DISTCORRALL    Calculate multivariate distance correlation between each pair of two random variable sets in X
% ==============================================================================================
% [ INPUTS ]
%     X = a 1D cell array of length m. Each cell contains 2D time-series matrices with n time points and p variates (voxels) (n x p)
%           p can be different across m, but n should be the same.
% .   optCent = Centering option
%                   1: U-centering (default) / 0: double-centering
% -----------------------------------------------------------------------------------------------
% [ OUTPUTS ]
%     dCor = m by m multivariate distance correlation matrix
%     dCov = m by m distance covariance
%     dVar = m distance variance
%     stats = independent t and p-value of distance correlation coefficient.
%
% Last update: Jan 31, 2019.
%
% Copyright 2019. Kwangsun Ray Yoo (K Yoo), PhD
%     E-mail: kwangsun.yoo@yale.edu / rayksyoo@gmail.com
%     Dept. of Psychology
%     Yale University
% ================================================================================================

if nargin < 2;    optCent = 1;   end;
m = length(X);    n = size(X{1}, 1);

for x = 1:m
    % Removing voxels of which values are zero all the time.  These zero value voxels are considered from outside the brain.
    if length(find(sum(X{x})==0))
        % nVox2remov(x,1) = length(find(sum(X{x})==0));
        X{x}(:,[find(sum(X{x})==0)]) = [];
    end
    
    % Euclidean distance among time points of nodes in X
    temp_Edist1d = pdist(X{x});

    % 1D distance array to 2D matrix
    n2 = ceil( sqrt(length(temp_Edist1d)*2) );
    temp_Edist = tril(ones( n2 ),-1);
    temp_Edist(temp_Edist==1) = temp_Edist1d;    clear temp_Edist1d
    temp_Edist = temp_Edist + temp_Edist';
    
    % Centering
    if optCent == 1
        % U-centering
        temp_EdistCent = temp_Edist - repmat(sum(temp_Edist,2)/(n-2), 1,n) - repmat(sum(temp_Edist,1)/(n-2),n,1) + repmat( sum(sum(temp_Edist))/((n-1)*(n-2)), n, n);
        temp_EdistCent(find(eye(n))) = 0;
    else
        % double-centering
        temp_EdistCent = temp_Edist - repmat(mean(temp_Edist,2), 1,n) - repmat(mean(temp_Edist,1),n,1) + repmat( mean(mean(temp_Edist)), n, n);
    end;
    Edist_centered(:,:,x) = temp_EdistCent;    clear temp_*
end

if optCent == 1
    K = n*(n-3);
elseif optCent == 2
    K = n^2;
end;

% Distance variance
dVar = squeeze(sum(sum(Edist_centered.^2,1),2))/K;

% Distance covariance
dCov = zeros(m);
for nd1 = 1:m-1
    for nd2 = nd1+1:m
        dCov(nd1,nd2) = sum(sum(Edist_centered(:,:,nd1) .* Edist_centered(:,:,nd2)))/K;
    end
end;
dCov = dCov + dCov';
dCov(find(eye(m))) = dVar;

% Distance correlation
dVarSqrt = sqrt(dVar * dVar');
dCor = sqrt( dCov ./ dVarSqrt);
dCor(dCov<=0) = 0;
dCor(logical(eye(size(dCor)))) = 0;

% Stats
stats.df = K/2 -1;
stats.T = sqrt(stats.df) * (dCor ./ sqrt(1-dCor.^2));
stats.p = tcdf(stats.T, stats.df, 'upper');
