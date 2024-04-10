function pt = findchangeptssimple(x)
n = length(x);
    for ii = 0:length(x)
        %Chooses a point and divides the signal into two sections.
        ind1 = 1:ii;
        ind2 = (ii+1):n;
        partition1 = x(ind1);
        partition2 = x(ind2);

        %Computes an empirical estimate of the desired statistical property for each section.
        % Linear fit
        beta1 = [ind1(:) , ones(ii,1)] \ partition1(:);
        beta2 = [ind2(:) , ones(n-ii,1)] \ partition2(:);

        %At each point within a section, measures how much the property deviates from the empirical estimate. Adds the deviations for all points.
        res1 = partition1(:) - [ind1(:) , ones(ii,1)] * beta1;
        res2 = partition2(:) - [ind2(:) , ones(n-ii,1)] * beta2;

        %Adds the deviations section-to-section to find the total residual error.
        totalres(ii+1) = sqrt(sum(res1.^2) + sum(res2.^2));
    end
    %Varies the location of the division point until the total residual error attains a minimum.
    [minresidual, minind] = min(totalres);
    pt = minind;
end