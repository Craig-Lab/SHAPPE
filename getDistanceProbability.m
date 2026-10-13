function distance = getDistanceProbability(distribution,cum_probs,prob)

    
    upperBounds = distribution(2:end);
    lowerBounds = distribution(1:end-1);
    bin = find(prob>cum_probs,1,'last');
    if(isempty(bin))
        if numel(distribution) >= 2
            bin = 1;
        else
            distance = distribution(1);
            return
        end
    end
    ll = lowerBounds(bin);
    rl = upperBounds(bin);
    distance = ll+(rl-ll)*rand();
        
end