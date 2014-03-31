function [C_with,C_without] = multiidx_box_set(shape,min_idx)

% [C_with,C_without] = MULTIIDX_BOX_SET(shape,min_idx)
%
% given an index shape, generates C_with, the box indeces set up to that i.
% min_idx is either 0 or 1. C_without is C_with without the shape_midx
% itself


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2014 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N=length(shape);
shape_fun=@(i) i-shape(1:length(i)); 
C_with=multiidx_gen(N,shape_fun,0,min_idx);

C_without=C_with;
C_without(end,:)=[];

