
% Jake Davis
% 10/07/2018

function [depFlag, newVerts, PSup] = supNewPanel(verts,P,cntr,Mach)

% Domain of Dependence Check
% check for fully out, fully in, or what is partially in
B = sqrt(Mach^2-1);
[depFlag, pntData] = supDoDCheck(B,verts,P,cntr);

% Build new panel

if depFlag == 0
% completely outside -> do nothing
else
    % Coord. transformations to panel local scaled system
    transMat = supCoordTrans(verts,Mach,0,0);
    PSup = (P - cntr)';
    % PSup = transMat * (P - cntr)';
    vertsSup(1,:) = transMat * (verts(1,:) - cntr)';
    vertsSup(2,:) = transMat * (verts(2,:) - cntr)';
    vertsSup(3,:) = transMat * (verts(3,:) - cntr)';
    % Plot scaled panel and Mach cone
    [~,~,cond] = triParams(vertsSup,PSup,1);
    supPlotMachCone(-1,PSup(1),100,PSup,Mach,1)
    xiMin = PSup(1) - sqrt(PSup(3)^2);
    plot(xiMin,PSup(2),'m*')
    title('Scaled CSYS')
    xlim([-3 3])
    ylim([-3 3])
    
    if depFlag == 1
        % completely inside
        newVerts = vertsSup;
    elseif depFlag == 2
        % pnt(s) inside -> find intersection(s)
        j = 1;
        if cond(2) == 1
        % (x,y) inside panel
            newVerts(j,:) = [PSup(1) - sqrt(PSup(3)^2), PSup(2) 0];
            j = j + 1;
        end
        loopVerts = vertsSup;
        loopVerts(end+1,:) = loopVerts(1,:);

        for i = 1:3
            % Check if ith vertex is inside
            if pntData(i) == 1
                newVerts(j,:) = vertsSup(i,:);
                j = j + 1;
            end
            % Check for edge-Mach cone intersections
            [pnt1,pnt2] = supIntersectionPnt(PSup,loopVerts(i,:),loopVerts(i+1,:));
            % Check if intersection point is on the edge
            if loopVerts(i,2) < loopVerts(i+1,2)
                if pnt1(2) > loopVerts(i,2) && pnt1(2) < loopVerts(i+1,2)
                    newVerts(j,:) = pnt1;
                    j = j + 1;
                end
                if pnt2(2) > loopVerts(i,2) && pnt2(2) < loopVerts(i+1,2)
                    newVerts(j,:) = pnt2;
                    j = j + 1;
                end
            else
                if pnt1(2) < loopVerts(i,2) && pnt1(2) > loopVerts(i+1,2)
                    newVerts(j,:) = pnt1;
                    j = j + 1;
                end
                if pnt2(2) < loopVerts(i,2) && pnt2(2) > loopVerts(i+1,2)
                    newVerts(j,:) = pnt2;
                    j = j + 1;
                end
            end
            hold all
            plot(pnt1(1),pnt1(2),'ob')
            plot(pnt2(1),pnt2(2),'ob')
        end
    elseif depFlag == 3
        % Check for edge-Mach cone intersections
        depFlag = 0;
        i = 1;
        loopVerts = vertsSup;
        loopVerts(end+1,:) = loopVerts(1,:);
        while i <= 3  && depFlag ~= 3
            [pnt1,pnt2] = supIntersectionPnt(PSup,loopVerts(i,:),loopVerts(i+1,:));
            if loopVerts(i,2) < loopVerts(i+1,2)
                if pnt1(2) > loopVerts(i,2) && pnt1(2) < loopVerts(i+1,2)
                    depFlag = 3;
                elseif pnt2(2) > loopVerts(i,2) && pnt2(2) < loopVerts(i+1,2)
                    depFlag = 3;
                end
            else
                if pnt1(2) < loopVerts(i,2) && pnt1(2) > loopVerts(i+1,2)
                    depFlag = 3;
                elseif pnt2(2) < loopVerts(i,2) && pnt2(2) > loopVerts(i+1,2)
                    depFlag = 3;
                end
            end
            i = i + 1;
        end
        newVerts = vertsSup;
    end
end

% Plot Stuff
plotVerts = newVerts;
plotVerts(end+1,:) = plotVerts(1,:);
figure
hold on
plot3(plotVerts(:,1), plotVerts(:,2), plotVerts(:,3),'-ok')
supPlotMachCone(-1,PSup(1),100,PSup,Mach,1)
title('Scaled Trimmed CSYS')
axis equal
grid on

end

%% Old intersection identification algorithm

% % Build new panel
% if depFlag == 0         % completely outside -> no influence
%     newVerts = vertsSup;
% elseif depFlag == 1     % completely inside -> do normal stuff
%     newVerts = vertsSup;
% elseif depFlag == 2     % pnt(s) inside -> find intersection(s)
%     if pntData(1) == 1 && pntData(2) == 1
%     % pnt1 & pnt2 are in -> edge2 & edge3 intersection
%     % edge2 -> new pnt3     % edge3 -> new pnt4
%         newVerts(3,:) = [supIntersectionPnt(PSup,vertsSup(2,:),vertsSup(3,:)), 0];
%         newVerts(4,:) = [supIntersectionPnt(PSup,vertsSup(3,:),vertsSup(1,:)), 0];
%     elseif pntData(2) == 1 && pntData(3) == 1
%     % pnt2 & pnt3 are in -> edge1 & edge3 intersection
%     % edge1 -> new pnt1     % edge3 -> new pnt4
%         newVerts(1,:) = [supIntersectionPnt(PSup,vertsSup(1,:),vertsSup(2,:)), 0];
%         newVerts(4,:) = [supIntersectionPnt(PSup,vertsSup(3,:),vertsSup(1,:)), 0];
%     elseif pntData(3) == 1 && pntData(1) == 1
%     % pnt3 & pnt1 are in -> edge1 & edge2 intersection
%     % edge1 -> new pnt2     % edge2 -> new pnt3
%         newVerts(2,:) = [supIntersectionPnt(PSup,vertsSup(1,:),vertsSup(2,:)), 0];
%         newVerts(4,:) = vertsSup(3,:);
%         newVerts(3,:) = [supIntersectionPnt(PSup,vertsSup(2,:),vertsSup(3,:)), 0];
%     elseif pntData(1) == 1
%     % pnt1 is in -> edge1 & edge3 intersection
%     % edge1 -> new pnt2     % edge3 -> new pnt3
%         newVerts(2,:) = [supIntersectionPnt(PSup,vertsSup(1,:),vertsSup(2,:)), 0];
%         newVerts(3,:) = [supIntersectionPnt(PSup,vertsSup(3,:),vertsSup(1,:)), 0];
%     elseif pntData(2) == 1
%     % pnt2 is in -> edge1 & edge2 intersection
%     % edge1 -> new pnt1     % edge2 -> new pnt3
%         newVerts(1,:) = [supIntersectionPnt(PSup,vertsSup(1,:),vertsSup(2,:)), 0];
%         newVerts(3,:) = [supIntersectionPnt(PSup,vertsSup(2,:),vertsSup(3,:)), 0];
%     elseif pntData(3) == 1
%     % pnt3 is in -> edge2 & edge3 intersection
%     % edge2 -> new pnt2     % edge3 -> new pnt1
%         newVerts(2,:) = [supIntersectionPnt(PSup,vertsSup(2,:),vertsSup(3,:)), 0];
%         newVerts(1,:) = [supIntersectionPnt(PSup,vertsSup(3,:),vertsSup(1,:)), 0];
%     end
% elseif depFlag == 3     % no pnt(s) inside -> don't need intersections
%     % do stuff to get w012 and whatever else
%     % or completely outside 
% end
