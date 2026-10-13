function s = sampleTargetDistance(Neigh,distributionsObs,distributionsRandom,copulaInformation,types,...
    samplingindexes,types_names,correlatedSamplingIdx,uncorrelatedSamplingIdx,exclusionDistance)
    
    copulaType = lower(string(copulaInformation.kernel));
    rho = copulaInformation.rho;
    if(copulaType == "t")
        nu = copulaInformation.nu;
    end

    Ntypes = length(types);
    s = zeros(1,Ntypes*Neigh);

    % sample a copula row
    if(copulaType=="t")
        U = copularnd('t',rho,nu,1);
    elseif(copulaType=="gaussian")
        U = copularnd('Gaussian',rho,1);
    else
        % No copula: sample indpendant probabilities
        U = unifrnd(0,1,[1 size(rho,2)]);
    end

    for nD=correlatedSamplingIdx
        % find which type this nD belongs to
        current_dist = [];
        for t = 1:Ntypes
            ct_name  = types_names(t);

            if nD >= samplingindexes(1,t) && nD <= samplingindexes(2,t)
                if(isfield(distributionsObs,ct_name))
                    current_dist = distributionsObs.(ct_name).distribution{nD - samplingindexes(1,t) + 1};
                    cumProb = distributionsObs.(ct_name).cum_probs{nD - samplingindexes(1,t) + 1};
    
                    if isempty(current_dist) || ~any(current_dist > 0)
                        current_dist = [];
                        break
                    end
    
                    if nD == samplingindexes(1,t)
                        cumulDistance = 0;
                    end
                end
                break
            end
        end

        if isempty(current_dist)
            continue
        end

        ind_U = find(correlatedSamplingIdx==nD, 1, "first");
        copulaInformation.sampleIndex = ind_U;

        foundDistance = false;
        while ~foundDistance
            if numel(correlatedSamplingIdx) > 1
                u = U(1, ind_U);
            else
                u = rand();
            end
            u = min(max(u, realmin), 1-realmin);
            d = getDistanceProbability(current_dist, cumProb,u);
            foundDistance = isfinite(d) && (d >= 0);
        end

        if nD < (Ntypes*Neigh + 1)
            cumulDistance = cumulDistance + d;
            d = cumulDistance;
        end

        s(1,nD) = d;
    end % End of for nD=correlatedSamplingIdx

    for nD=uncorrelatedSamplingIdx
        current_dist = [];

        for t = 1:Ntypes
            ct_name  = types_names(t);

            if nD >= samplingindexes(1,t) && nD <= samplingindexes(2,t)
                current_dist = distributionsRandom.(ct_name).distribution{nD - samplingindexes(1,t) + 1};
                cumProb = distributionsRandom.(ct_name).cum_probs{nD - samplingindexes(1,t) + 1};

                if isempty(current_dist) || ~any(current_dist > 0)
                    current_dist = [];
                    break
                end

                if nD == samplingindexes(1,t)
                    cumulDistance = 0;
                end
                break
            end
        end

        if isempty(current_dist)
            continue
        end

        foundDistance = false;
        while ~foundDistance
            u = rand();
            u = min(max(u, realmin), 1-realmin);
            d = getDistanceProbability(current_dist, cumProb,u);
            foundDistance = isfinite(d) && (d >= 0);
        end

        cumulDistance = cumulDistance + d;
        s(1,nD) = cumulDistance;
    end % End of for nD=uncorrelatedSamplingIdx

end