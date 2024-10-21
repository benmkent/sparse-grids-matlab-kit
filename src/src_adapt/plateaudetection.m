function [prof_tol, changeptind, beta0,beta1, detect_flag] = plateaudetection(inds, env, controls)
% PLATEAUDETECTION Detects a plateau in the input map inds->log10(env)
%
% This is used to detect the presence of noise in a spectral representation
% of a function for polynomial total degree represented by inds, and
% coefficient norms envelope env.
% This assumes that the logarithm of the spectral envelope can be segmented into two linear models (representing a decay in spectral coefficients and the spectral coefficient plateau)
%
% [prof_tol, changeptind, beta0,beta1, detect_flag] = PLATEAUDETECTION(inds, env, controls)
% Inputs:
%   inds: index vector for envelope
%   env: envelope to detect plateau in
%   controls: plateau detection controls
%
% Outputs:
%   prof_tol:       detected plateau level
%   changptind:     index of start of plateau
%   beta0:          model beta0[1] + beta0[2]*i for spectral decay
%   beta1:          model beta1[1] + beta1[2]*i for spectral plateau
%   detect_flag:    flag indicating a plateau is identified.
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

%%  Trim envelope
% for an analytic function we expect coefficient norms to decrease
% at an exponential rate with respect to the total degree of the
% polynomials.
% This may not be initially so we use a control burn_in.
% The final coefficient is often not representative. This can be
% ignored with the burn out

offset = controls.burn_in;
tail = controls.burn_out;

% change point detection is used to partition log10(envelope) into
% two --- hopefully initially a linear (i.e. exponential envelope)
% region followed by a constant plateau.
log10env = log10(env((1+offset):(end-tail))+eps()); % use eps to avoid issues if envelope goes to zero.

%% Detect change point
% Provided trimmed envelope is suitably long, apply change point
% detection to the log10env.
% We seek a splitting into two regiemes.
min_length = controls.min_plateau_length;
if length(log10env) >= min_length
    % Change point detection on log10env
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave || ~exist('findchangepts')
        ipt = findchangeptssimple(log10env);
    else
        ipt = findchangepts(log10env,Statistic="linear",MaxNumChanges=1);
    end
else
    ipt = [];
end

%% Compute the linear functions
% compute first linear function
inds_0 = inds(offset+1:(offset+ipt-1));
beta0 = [ones(length(inds_0),1), inds_0.'] \ log10env(1:(ipt-1)).';
% compute the linear function approximating the second region
inds_1 = inds(offset+ipt:end-tail);
beta1 = [ones(length(inds_1),1), inds_1.'] \ log10env(ipt:(end)).';

length_plateau = length(log10env) - ipt;

%% determine if the second region is a plateau
% here we test the gradient against a tolerance see if it is "flat"
% and if it is long enough
if ~all(beta1==0) && abs(beta1(2)) < controls.plateau_gradient ...
        && length_plateau >= min_length

    % If plateau the profit tolerance is set equal to the start of
    % the elbow.
    prof_tol = 10^(beta1(1) +inds(ipt)*beta1(2));
    prof_tol = prof_tol * controls.safety_factor;
    detect_flag = 1;
else
    % Else it is is reassigned to the initial tolerance.
    prof_tol = controls.tol_init;
    detect_flag = 0;
end

%% Correct ipt for offset and as index
changeptind = inds(ipt+offset);
