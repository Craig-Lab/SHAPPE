function symmetricsPoints = getSymmetricPoints(conf,R,maxDistance,shape)
    symmetricsPoints = conf;
    idx_symmetric_points = 1;
    if(strcmp(shape,"circle"))
        % Loop through points
        for i=1:size(conf,1)
            distance_from_border = R-norm(conf(i,1:2));
            if(distance_from_border<=maxDistance)
                rc = conf(i,1:2)./norm(conf(i,1:2));
                symmetricsPoints(idx_symmetric_points,:) = [conf(i,1:2)-2.*R.*rc,conf(i,3:end)];
                idx_symmetric_points = idx_symmetric_points+1;
            end
        end
    elseif(strcmp(shape,"square"))
        maxDistance = sqrt(2.*(maxDistance).^2);
        symmetricsPoints = [symmetricsPoints;symmetricsPoints;symmetricsPoints];
        y_up = R;
        y_down = -R;
        x_left = -R;
        x_right = R;
        displacements = [-1 0;0 -1;-1 -1];
        % Loop through points
        for i=1:size(conf,1)
            x = conf(i,1);
            y = conf(i,2);
            displacements = [-1 0;0 -1;-1 -1];
            if((y_up-y)<=maxDistance || (y-y_down)<=maxDistance || ...
                    (x-x_left)<=maxDistance || (x_right-x)<=maxDistance)
                if(x<0)
                    displacements(:,1) = -1.*displacements(:,1);
                end
                if(y<0)
                    displacements(:,2) = -1.*displacements(:,2);
                end
                for d=1:3
                    new_x = conf(i,1)+2.*R.*displacements(d,1);
                    new_y = conf(i,2)+2.*R.*displacements(d,2);
                    if(new_x>=(x_left-maxDistance) && new_x <= (x_right+maxDistance) && ...
                        new_y>=(y_down-maxDistance) && new_y<=(y_up+maxDistance))
                        symmetricsPoints(idx_symmetric_points+d-1,:) = [new_x,new_y,conf(i,3:end)];
                    end
                end
                idx_symmetric_points = idx_symmetric_points+3;
            end
        end
    end
    symmetricsPoints(idx_symmetric_points:end,:) = [];
end