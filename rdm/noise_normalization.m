function segs=noise_normalization(segs,varargin)
% Perform noise normalization before calculating the RDM
% Based on recommendations from: 
%       Guggenmos, Sterzer, and Cichy. "Multivariate pattern analysis for MEG:
%       A comparison of dissimilarity measures." NeuroImage (2018).
%
% Default settings: multivariate normalization with covariance calculated
%   based on the full segment, and no separation by categories.
%
% Example usage: segs_normalized =
%   noise_normalization(segs,'univariate',1:25,'categories',cat_vec)
% Performs noise normalization on the segments, with diagonal covariance
% (univariate) based on time-points 1:25 of the segment and separately for
% each category based on the numbers in cat_vec.
%
% Input:
%   segs (segments for RSA, dimensions: n_stim x n_features x n_time)
%   Optional:
%       - 'multivariate'\'univariate' - method for calculating the
%           covariance matrix. multivariate = full covariance (default),
%           univariate = only diagonal (single feature variances).
%       - 'per_point' - calculate the covariance per time-point, i.e. each 
%           time-point is normalized separately (slower) (overrides
%           specific indices if they are given)
%       - indices_vector (numeric vector) - time-points to calculate the 
%           covariance matrix (indices in segs). Not relevant in case 'per_point'.
%           Default is the entire segment.
%       - In case the data comes from several categories use: 'categories'
%           followed by numeric vector with values between 1:n_cat to indicate
%           which trial belongs to which category. Covariances will be calculated
%           per category then averaged according to the respective sizes.
%
% Output: segs (same dimensions), after noise normalization
%
% Multivariate covariance is calculated using the Ledoit-Wold formula (function cov1para.m)
% Ledoit, and Wolf. "A well conditioned estimator for large dimensional covariance matrices." (2000).
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2020
% Send bug reports and requests to gal.vishne@gmail.com

n_stim = size(segs,1); n_features = size(segs,2); n_time = size(segs,3);

multivariate = true;
per_point = false;
norm_range = 1:n_time;    %default is full epoch
categories = ones(1,n_stim);  %default is one category
arg  = 1;
while arg <= size(varargin,2)
    if isnumeric(varargin{arg})
        norm_range = varargin{arg};
        if norm_range(1) < 1 || norm_range(end) > n_time
            error('Given range exceeds segment duration');
        end
        if size(norm_range,1)~=1
            norm_range = norm_range';
        end
        arg = arg + 1;
    else
        switch varargin{arg}
            case 'multivariate'
                multivariate = true;
                arg = arg + 1;
            case 'univariate'
                multivariate = false;
                arg = arg + 1;
            case 'per_point'
                per_point = true;
                arg = arg + 1;
            case 'categories'
                categories = varargin{arg+1};
                if length(categories) ~= size(segs,1)
                    error('Category vector is not the same length as the number of segments')
                end
                arg = arg + 2;
            otherwise
                error(['Unknown optional argument name: ' varargin{arg} '.']);
        end
    end
end

% don't weight tiny categories the same as the big ones
category_weights = nan(1,max(categories));
for c=1:max(categories)
    category_weights(c) = (sum(categories==c)/length(categories))*max(categories);
end

if ~per_point
    sigma = nan(n_features,n_features,length(norm_range), max(categories));
    for c = 1:max(categories)
        for t = norm_range
            if multivariate
                sigma(:,:,norm_range==t,c) = category_weights(c)*cov1para(segs(categories==c,:,t));
            else
                sigma(:,:,norm_range==t,c) = category_weights(c)*diag(diag(cov(segs(categories==c,:,t))));
            end
        end
    end
    sigma = mean(sigma,[3 4])^-0.5;

    for t = 1:n_time
        segs(:,:,t) = segs(:,:,t) * sigma;
    end
else
    for t = 1:n_time
        sigma = nan(n_features,n_features,max(categories));
        for c=1:max(categories)
            if multivariate
                sigma(:,:,c) = category_weights(c)*cov1para(segs(categories==c,:,t));
            else
                sigma(:,:,c) = category_weights(c)*diag(diag(cov(segs(categories==c,:,t))));
            end
        end
        sigma = mean(sigma,3)^-0.5;
        
        segs(:,:,t) = segs(:,:,t) * sigma;
    end
end
    
end