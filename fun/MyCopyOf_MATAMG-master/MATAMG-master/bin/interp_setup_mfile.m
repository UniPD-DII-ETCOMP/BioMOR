function [maxindex carray] = interp_setup_mfile(C,ia,carray)
maxindex = 1;

m = length(ia)-1;
% Measure number of nonzeros of P and define map from coarse to fine grid
for i = 1 : m
    % at row corresponding to C point, P has only one nonzero.
    if C(i) 
        maxindex = maxindex + 1;
    % maximum possible number of nonzero at F point is its nbd point.
    else
        maxindex = maxindex + ia(i+1) - ia(i);
    end
    
    % Find coarse level index for C points
    if i == 1
        carray(i) = C(i);
    else
        carray(i) = carray(i-1) + C(i);
    end
end


