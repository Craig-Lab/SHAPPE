% Place points in points_index with an exclusion distance
function config = placeWithExclusion(config, points_index, placedTypes, Regions, options, R,maxTries)

    % Find points that are already placed.
    placed_points = find(ismember(config(:,end), placedTypes));
    placed_points_positions = config(ismember(config(:,end),placedTypes),1:2);
    for nC=(points_index)'
        pos = config(nC,1:2);
        evaluatePosition = true;
        numTrials = 0;
        best_pos = pos;
        current_min_distance = -inf;
        while(evaluatePosition)
            numTrials = numTrials+1;
            if isempty(placed_points)
                ok = true;
                mind2 = inf;
            else
                d2 = sum((placed_points_positions - pos).^2, 2);
                ok = all(d2 >= options.exclusionDistance.^2);
                mind2 = min(d2);
            end
    
            if(~options.useRegions)
                okRegions = true;
            else
                if(options.attract)
                    okRegions = false;
                    if(~isempty(Regions))
                        for vv=1:size(Regions,1)
                            dv = norm(pos-[Regions(vv,1),Regions(vv,2)]);
                            if(dv<Regions(vv,3))
                                okRegions = true;
                                break
                            end % End of if(dv<Regions(vv,3))
                        end % End of for vv=1:size(Regions,1)
                    else
                        okRegions = true;
                    end % End of if(~isempty(Regions) && numTrials>1) ... else ... end
                else
                    okRegions = true;
                    if(~isempty(Regions))
                        for vv=1:size(Regions,1)
                            dv = norm(pos-[Regions(vv,1),Regions(vv,2)]);
                            if(dv<Regions(vv,3))
                                okRegions = false;
                                break
                            end % End of if(dv<Regions(vv,3))
                        end % End of for vv=1:size(Regions,1)
                    end % End of if(~isempty(Regions) && numTrials>1) ... else ... end
                end
            end
            if(~ok && mind2>current_min_distance && okRegions)
                current_min_distance = mind2;
                best_pos = pos;
            end
            if(ok && okRegions)
                evaluatePosition = false;
                config(nC,1:2) = pos;
                placed_points(end+1,1) = nC;
                placed_points_positions(end+1,:) = pos;
            elseif(numTrials>=maxTries)
                evaluatePosition = false;
                pos = best_pos;
                config(nC,1:2) = pos;
                placed_points(end+1,1) = nC;
                placed_points_positions(end+1,:) = pos;
            else
                % Generate new random position
                if(strcmp(options.shapeSynthetic,"square"))
                    pos = [unifrnd(-R,R) unifrnd(-R,R)];
                else
                    pos_theta = 2.*pi.*rand();
                    pos_r = R*sqrt(rand());
                    pos_x = pos_r.*cos(pos_theta);
                    pos_y = pos_r.*sin(pos_theta);
                    pos = [pos_x pos_y];
                end
            end
        end
    end 

end