function [libDist,shape] = interpProfile5(imgA,imgB,divide)

% to extract the point from image.
%[pA, pB] = extractPoint(imgA,imgB);
 pA = imgA;
 pB = imgB;
% % to reduce the set of points from image.
%try
%set(0,'RecursionLimit',600);
%[psAsim,iax] = dpsimplify(pA,1);
%[psBsim,ibx] = dpsimplify(pB,1);


% catch
%     psAsim = imgA;
%     psBsim = imgB;
%     iax = 1:1:size(psAsim,2);
%     ibx = 1:1:size(psBsim,2);
% end
psA = pA;
psB = pB;

% sA = size(psA,1); %point in A
% sB = size(psB,1); %point in B
sAsim = size(psA,1); %point in A
sBsim = size(psB,1); %point in B
if sAsim < sBsim
    
    temp = psA;
    psA = psB;
    psB = temp;
    
%     temp = sA;
%     sA = sB;
%     sB = temp;
    
end

figure,plot(psA(:,1),psA(:,2),'or');
hold on,plot(psB(:,1),psB(:,2),'xb');
%hold on,plot(psAsim(:,1),psAsim(:,2),'.g');
%hold on,plot(psBsim(:,1),psBsim(:,2),'.y');
hold off;

% psA = psAsim;
% psB = psBsim;
sA = size(psA,1); %point in A
sB = size(psB,1); %point in B

matching = knnclassify(psA, psB, [1:1:sB]');
psBtemp = [];
for teemp = sB:-1:1
    if sB == matching
        continue;
    else
        psBtemp =[psBtemp; psB(teemp,1:2)];
        
    end
end
sBt = size(psBtemp,1);
matching2 = knnclassify(psA,psBtemp, [1:1:sBt]');
%disp(['Hello']);
new_out = zeros(sA,2);
divide = abs(divide);

libDist = [];
%divide = 30;
for step = 1:divide
    for pi = 1:1:sA
        try
            if (new_out(pi,1)== psB(matching(pi),1)) && (new_out(pi,2) == psB(matching(pi),2))
                continue;   % decrease loop for corresponsedense point
                            % made lastest profile again
            end
            
          
            distX = pdist([psA(pi,1) ; psB(matching(pi),1)],'euclidean');
            distY = pdist([psA(pi,2) ; psB(matching(pi),2)],'euclidean');
            %% x and y is equal
            if (psA(pi,1) == psB(matching(pi),1)) && (psA(pi,2) == psB(matching(pi),2))
                new_out(pi,1) = psA(pi,1) ;
                new_out(pi,2) = psA(pi,2) ;
                
            %% x is different
            elseif psA(pi,1) ~= psB(matching(pi),1)
                if psA(pi,1) > psB(matching(pi),1)
                    nX =psA(pi,1)-(step*(distX/divide));
                    if nX < psB(matching(pi),1)
                        nX = psB(matching(pi),1);
                    end
                elseif psA(pi,1) < psB(matching(pi),1)
                    nX =psA(pi,1)+(step*(distX/divide));
                    if nX > psB(matching(pi),1)
                        nX = psB(matching(pi),1);
                    end
                end
                nY = interp1([psA(pi,1);psB(matching(pi),1)],[psA(pi,2);psB(matching(pi),2)],nX,'spline');
                new_out(pi,1) = nX;
                new_out(pi,2) = nY;
            %% x is not different
            else
                if psA(pi,2) > psB(matching(pi),2)
                    nX =psA(pi,2)-(step*(distY/divide));
                    if nX < psB(matching(pi),2)
                        nX = psB(matching(pi),2);
                    end
                elseif psA(pi,2) < psB(matching(pi),2)
                    nX =psA(pi,2)+(step*(distY/divide));
                    if nX > psB(matching(pi),2)
                        nX = psB(matching(pi),2);
                    end
                end
                nY = interp1([psA(pi,2);psB(matching(pi),2)],[psA(pi,1);psB(matching(pi),1)],nX,'spline');
                new_out(pi,1) = nY;
                new_out(pi,2) = nX;
            end
            
            if isnan(nY)
               disp(['ffff']) 
            end
        catch
        end
        % =============== width condition ==================
        if pi == 1
            P1 = new_out(pi,1:2);
        elseif pi == ceil((sA-1)/2) % ======= experiment ======= 
            % half of shape 
            P2 = new_out(pi,1:2);
            D = pdist([P1;P2],'euclidean');
        end   
    end
    
    libDist = [libDist;D];
    if step == divide
    %close all
        figure,plot3(psA(:,1),psA(:,2),ones(sA,1),'r*-');
        hold on,plot3(psB(:,1),psB(:,2),ones(sB,1)*2,'b*-');
        hold on,plot3(new_out(:,1),new_out(:,2),ones(sA,1)*3,'g*-');

        title(['step : ' num2str(step) ' size A = ' num2str(sA) ' size out = ' num2str(size(new_out,1))])
    end
    %pause(0.5)
    
%     recon1 = repmat(new_out(:,1)',5,1);
%     recon2 = repmat(new_out(:,2)',5,1);
%     recon3 = 1:5;
%     recon3 = recon3';
%     recon3 = repmat(recon3,1,48);
%     figure,surfnorm(recon1,recon2,recon3);
    shape(step).XY = new_out;
end
%%
new_out2 = zeros(sA,2);
libDist2 = [];
for step = 1:divide
    for pi = 1:1:sA
        try
            if (new_out2(pi,1)== psBtemp(matching2(pi),1)) && (new_out2(pi,2) == psBtemp(matching2(pi),2))
                continue;   % decrease loop for corresponsedense point
                            % made lastest profile again
            end       
            distX = pdist([psA(pi,1) ; psBtemp(matching2(pi),1)],'euclidean');
            distY = pdist([psA(pi,2) ; psBtemp(matching2(pi),2)],'euclidean');
            %% x and y is equal
            if (psA(pi,1) == psBtemp(matching2(pi),1)) && (psA(pi,2) == psBtemp(matching2(pi),2))
                new_out2(pi,1) = psA(pi,1) ;
                new_out2(pi,2) = psA(pi,2) ;
                
            %% x is different
            elseif psA(pi,1) ~= psBtemp(matching2(pi),1)
                if psA(pi,1) > psBtemp(matching2(pi),1)
                    nX =psA(pi,1)-(step*(distX/divide));
                    if nX < psBtemp(matching2(pi),1)
                        nX = psBtemp(matching2(pi),1);
                    end
                elseif psA(pi,1) < psBtemp(matching2(pi),1)
                    nX =psA(pi,1)+(step*(distX/divide));
                    if nX > psBtemp(matching2(pi),1)
                        nX = psBtemp(matching2(pi),1);
                    end
                end
                nY = interp1([psA(pi,1);psBtemp(matching2(pi),1)],[psA(pi,2);psBtemp(matching2(pi),2)],nX,'spline');
                new_out2(pi,1) = nX;
                new_out2(pi,2) = nY;
            %% x is not different
            else
                if psA(pi,2) > psBtemp(matching2(pi),2)
                    nX =psA(pi,2)-(step*(distY/divide));
                    if nX < psBtemp(matching2(pi),2)
                        nX = psBtemp(matching2(pi),2);
                    end
                elseif psA(pi,2) < psBtemp(matching2(pi),2)
                    nX =psA(pi,2)+(step*(distY/divide));
                    if nX > psBtemp(matching2(pi),2)
                        nX = psBtemp(matching2(pi),2);
                    end
                end
                nY = interp1([psA(pi,2);psBtemp(matching2(pi),2)],[psA(pi,1);psBtemp(matching2(pi),1)],nX,'spline');
                new_out2(pi,1) = nY;
                new_out2(pi,2) = nX;
            end
            
            if isnan(nY)
               disp(['ffff']) 
            end
        catch
        end
        % =============== width condition ==================
        if pi == 1
            P1 = new_out2(pi,1:2);
        elseif pi == ceil((sA-1)/2) % ======= experiment ======= 
            % half of shape 
            P2 = new_out2(pi,1:2);
            D = pdist([P1;P2],'euclidean');
        end   
    end
    
    libDist2 = [libDist;D];
    if step == divide
        %close all
        figure,plot3(psA(:,1),psA(:,2),ones(sA,1),'r*-');
        hold on,plot3(psBtemp(:,1),psBtemp(:,2),ones(sBt,1)*2,'b*-');
        hold on,plot3(new_out2(:,1),new_out2(:,2),ones(sA,1)*3,'g*-');

        title(['step : ' num2str(step) ' size A = ' num2str(sA) ' size out = ' num2str(size(new_out2,1))])
        %pause(0.5)
    end
%     recon1 = repmat(new_out(:,1)',5,1);
%     recon2 = repmat(new_out(:,2)',5,1);
%     recon3 = 1:5;
%     recon3 = recon3';
%     recon3 = repmat(recon3,1,48);
%     figure,surfnorm(recon1,recon2,recon3);
    shape1(step).XY = new_out2;
end
%%
% preparing model
modelling(shape1);
end

function [pAnew, pBnew] = extractPoint(imgA,imgB)

% [pointAX pointAY] = find(imgA == 1);
% [pointBX pointBY] = find(imgB == 1);
 [pointAX] = imgA(:,1);[pointAY] = imgA(:,2);
 [pointBX] = imgB(:,1);[pointBY] = imgB(:,2);
[ pAnew ] = sortedPoint( [pointAX pointAY] );
[ pBnew ] = sortedPoint( [pointBX pointBY] );

end

function [ pout ] = sortedPoint( p )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

px = p(:,1);
py = p(:,2);

cx = mean(px);
cy = mean(py);

a = atan2(py-cy,px-cx);

[~,ord] = sort(a);

px = px(ord);
py = py(ord);

pout = [px py];

end
function [] = modelling(shape)
reconX = [];
reconY = [];
reconZ = [];
for j = 1:size(shape,2)
    recon = shape(j).XY;
    reconX = [reconX ;recon(:,1)'];
    reconY = [reconY ;recon(:,2)'];
    reconZ = [reconZ ;(ones(1,size(reconX,2))*j)];
    
end
%creating model
figure,hSurface = surf(reconX,reconY,reconZ,'EdgeColor','k','MeshStyle','column','FaceColor','interp','FaceLighting','gouraud');
set(hSurface,'FaceColor','c');
hold on;
plot3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','k-');
%fill3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','c')
plot3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','k-');
%fill3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','c');
axis tight;
end
