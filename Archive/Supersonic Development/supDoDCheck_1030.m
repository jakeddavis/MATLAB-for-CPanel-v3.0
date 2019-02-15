
% Jake Davis
% 10/14/2018

% Pretty sure all that matters is that the Mach cone intersects the panel
% in some way. If it does, move ahead with inf. coeff. calcs, and check for
% special cases. The special cases should catch what is causing atan2 to
% break

% Dependence Flag info
% 0 = completely outside
% 1 = completely inside
% 2 = point(s) inside
% 3 = pure edges

% Edge Data Info
% Ex: [1 0 1] means edges 1 and 3 are intersected, 2 is outside


function [depFlag, edgeData] = supDoDCheck_1030(M,verts,P,cntr,a,b)

B = sqrt(M^2-1);
pntInterData = zeros(3,1);
edgeData = zeros(3,1);

% c_r = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b)]';

% disp('P before')
% disp(P)
% disp('verts before')
% disp(verts)
% transWind = supCoordTrans(verts,M,a,b);
% % transWind = global2wind(a,b);
% P = (transWind * (P-cntr)')';
% verts(1,:) = transWind * (verts(1,:)-cntr)';
% verts(2,:) = transWind * (verts(2,:)-cntr)';
% verts(3,:) = transWind * (verts(3,:)-cntr)';
% disp('P after')
% disp(P)
% disp('verts after')
% disp(verts)
% [~,~,~] = triParams(verts,P,1);
% % supPlotMachConeZ(-1,P(1),100,P,verts,M,0)
% supPlotMachCone(-1,P(1),100,P,M,0)

% edgeFlag = supEdgeCheckZ(verts(1,:),verts(2,:),P,B)
% edgeFlag = supEdgeCheckZ(verts(2,:),verts(3,:),P,B)
% edgeFlag = supEdgeCheckZ(verts(3,:),verts(1,:),P,B)

% disp('P before')
% disp(P)
% disp('verts before')
% disp(verts)
% transLocal = getLocalSys(verts,1);
% P = (transLocal * (P-cntr)')';
% verts(1,:) = transLocal * (verts(1,:)-cntr)';
% verts(2,:) = transLocal * (verts(2,:)-cntr)';
% verts(3,:) = transLocal * (verts(3,:)-cntr)';
% disp('P after')
% disp(P)
% disp('verts after')
% disp(verts)
% [~,~,~] = triParams(verts,P,1);
% % supPlotMachConeZ(-1,P(1),100,P,verts,M,0)
% supPlotMachCone(-1,P(1),100,P,M,0)



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
                edgeData(2) = supEdgeCheck(verts(2,:),verts(3,:),P,B);
            elseif pntInterData(2) == 1
                edgeData(1) = 1;
                edgeData(2) = 1;
                edgeData(3) = supEdgeCheck(verts(3,:),verts(1,:),P,B);
%                 edgeData(3) = supEdgeCheckZ(verts(3,:),verts(1,:),P,B);
            elseif pntInterData(3) == 1
                edgeData(2) = 1;
                edgeData(3) = 1;
                edgeData(1) = supEdgeCheck(verts(1,:),verts(2,:),P,B);
            end
        elseif count == 2
            depFlag = 2;
            edgeData(:,:) = 1;
        else
            depFlag = 3;
            loopVerts = verts;
            loopVerts(end+1,:) = loopVerts(1,:);
            edgeDataZeros = zeros(3,1);
            for i = 1:3
                edgeData(i) = supEdgeCheck(loopVerts(i,:),loopVerts(i+1,:),P,B);
%                 edgeData(i) = supEdgeCheckZ(loopVerts(i,:),loopVerts(i+1,:),P,B);
            end
            
            if edgeData == edgeDataZeros
                depFlag = 0;
            end
        end
    end
end

end

%%

%             if pntInterData(1) == 1 && pntInterData(2) == 1
%                 edgeData(2) = 1;
%                 edgeData(3) = 1;
%             elseif pntInterData(2) == 1 && pntInterData(3) == 1
%                 edgeData(1) = 1;
%                 edgeData(3) = 1;
%             elseif pntInterData(1) == 1 && pntInterData(3) == 1
%                 edgeData(1) = 1;
%                 edgeData(2) = 1;
%             end

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
