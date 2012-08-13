function [C_without,C_with] = multiindex_box_set(shape,min_idx)

% [C_without,C_with] = MULTIINDEX_BOX_SET(shape,min_idx)
%
% given an index shape, generates C_with, the box indeces set up to that i.
% min_idx is either 0 or 1

N=length(shape);
shape_fun=@(i) i-shape(1:length(i)); 
C_with=multiidx_gen(N,shape_fun,0,min_idx);

C_without=C_with;
C_without(end,:)=[];

