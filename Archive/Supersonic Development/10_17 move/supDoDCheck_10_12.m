
% Jake Davis
% 10/09/2018

% Pretty sure all that matters is that the Mach cone intersects the panel
% in some way. If it does, move ahead with inf. coeff. calcs, and check for
% special cases. The special cases should catch what is causing atan2 to
% break

% Dependence Flag info
% 0 = completely outside -> no influence
% 1 = completely inside

function [depFlag, edgeData] = supDoDCheck_10_12(M,verts,P,cntr)

B = sqrt(M^2-1);
pntInterData = zeros(3,1);
edgeData = zeros(3,1);

% Define Stuff
depFlag = 0;
x = P(1); y = P(2); z = P(3);
x0 = cntr(1); y0 = cntr(2); z0 = cntr(3);
RBcntr = (x-x0)^2 - B^2*(y-y0)^2 - B^2*(z-z0)^2;
c0 = [1 0 0];

% compute stuff
Pc = P - P;     % set origin at point of interest (POI)
vertsCheck = verts - P;     % adjust panel to new CSYS

[~,Q,cond] = triParams(vertsCheck,Pc,0);
x0 = dot(Q-Pc,c0);
y0 = sqrt(norm(Q-Pc)^2 - x0^2);
if y0 >= B * x0
    dist = sqrt((x0+B*y0)^2/(1+B^2));
else
    dist = sqrt(norm(Q - Pc)^2);
end

if dot(P - cntr, c0) >= 0 || dist < cond(1)
    if dist > cond(1)
    % cntr of panel is further from domain than panel rad.
        if sign(RBcntr) == 1
        % cntr is completely inside domain
            depFlag = 1;
            edgeData(:,:) = 1;
        %%%%%%%%%%%%% don't need
        elseif sign(RBcntr) == -1
        % cntr is outside domain
        else
            error('uh oh...')
        %%%%%%%%%%%%% don't need
        end
    else
        count = 0;
        for i = 1:3
            xi = verts(i,1);
            yi = verts(i,2);
            zi = verts(i,3);
            RBverts = (x-xi)^2 - B^2*(y-yi)^2 - B^2*(z-zi)^2;
            if sign(RBverts) == 1 && dot(P - [xi yi zi], c0) >= 0
                count = count + 1;
                pntInterData(i) = 1;
            else
                pntInterData(i) = 0;
            end
        end
        
        if count == 3
            depFlag = 1;
            edgeData(:,:) = 1;
        elseif count == 1
            depFlag = 2;
            if pntInterData(1) == 1
                edgeData(1) = 1;
                edgeData(3) = 1;
            elseif pntInterData(2) == 1
                edgeData(1) = 1;
                edgeData(2) = 1;
            elseif pntInterData(3) == 1
                edgeData(1) = 2;
                edgeData(2) = 3;
            end
        elseif count == 2
            depFlag = 2;
            if pntInterData(1) == 1 && pntInterData(2) == 1
                edgeData(2) = 1;
                edgeData(3) = 1;
            elseif pntInterData(2) == 1 && pntInterData(3) == 1
                edgeData(1) = 1;
                edgeData(3) = 1;
            elseif pntInterData(1) == 1 && pntInterData(3) == 1
                edgeData(1) = 1;
                edgeData(2) = 1;
            end
        else
            depFlag = 3;
%             disp('pure edge intersection, or outside')
        end
    end
    
    if depFlag == 3
        % Check for edge-Mach cone intersections
        transMat = supCoordTrans(verts,M,0,0);
        PSup = transMat * (P - cntr)';
        supVerts(1,:) = transMat * (verts(1,:) - cntr)';
        supVerts(2,:) = transMat * (verts(2,:) - cntr)';
        supVerts(3,:) = transMat * (verts(3,:) - cntr)';
        
        [~,~,~] = triParams(supVerts,PSup,1);
        supPlotMachCone(-1,PSup(1),100,PSup,M,1)
        
        loopVerts = supVerts;
        loopVerts(end+1,:) = loopVerts(1,:);
        edgeInterDataZeros = zeros(3,1);
        for i = 1:3
            [pnt1,pnt2] = supIntersectionPnt(PSup,loopVerts(i,:),loopVerts(i+1,:));
            dot(P - pnt1, c0)
            dot(P - pnt2, c0)

            if dot(P - pnt1, c0) >= 0
                if loopVerts(i,2) < loopVerts(i+1,2)
                    if pnt1(2) > loopVerts(i,2) && pnt1(2) < loopVerts(i+1,2)
                        edgeData(i) = 1;
                    end
                else
                    if pnt1(2) < loopVerts(i,2) && pnt1(2) > loopVerts(i+1,2)
                        edgeData(i) = 1;
                    end
                end
            elseif dot(P - pnt2, c0) >= 0
                if loopVerts(i,2) < loopVerts(i+1,2)
                    if pnt2(2) > loopVerts(i,2) && pnt2(2) < loopVerts(i+1,2)
                        edgeData(i) = 1;
                    end
                else
                    if pnt2(2) < loopVerts(i,2) && pnt2(2) > loopVerts(i+1,2)
                        edgeData(i) = 1;
                    end
                end
            end
        end
        if edgeData == edgeInterDataZeros
            depFlag = 0;
        end
    end
end

end

%%

%         i = 1;
%         loopVerts = vertsSup;
%         loopVerts(end+1,:) = loopVerts(1,:);
%         while i <= 3  && depFlag == 0
%             [pnt1,pnt2] = supIntersectionPnt(P,loopVerts(i,:),loopVerts(i+1,:));
%             if loopVerts(i,2) < loopVerts(i+1,2)
%                 if pnt1(2) > loopVerts(i,2) && pnt1(2) < loopVerts(i+1,2)
%                     depFlag = 1;
%                 elseif pnt2(2) > loopVerts(i,2) && pnt2(2) < loopVerts(i+1,2)
%                     depFlag = 1;
%                 end
%             else
%                 if pnt1(2) < loopVerts(i,2) && pnt1(2) > loopVerts(i+1,2)
%                     depFlag = 1;
%                 elseif pnt2(2) < loopVerts(i,2) && pnt2(2) > loopVerts(i+1,2)
%                     depFlag = 1;
%                 end
%             end
%             i = i + 1;
%         end

%%

%         i = 1;
%         while i <=3 && depFlag == 0
%             xi = verts(i,1);
%             yi = verts(i,2);
%             zi = verts(i,3);
%             RBverts = (x-xi)^2 - B^2*(y-yi)^2 - B^2*(z-zi)^2;
%             if sign(RBverts) == 1
%             % panel corner is in domain
%                 depFlag = 2;
%             end
%             i = i + 1;
%         end
