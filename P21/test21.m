function [out] = test21()
%clear all;close all;
%% interpolation process
load('dataAlltoothpasteTest.mat');

%% rotation module
% temp = xy1Sh(:,1);
% xy1Sh(:,1) = xy1Sh(:,2);
% xy1Sh(:,2) = temp(:);

[libDist, shape] = interpProfile4(xy1Sh,xy2Sh,divide);
hold on;
%% sweep process

% Library pickup
model(1).XY = shape(1).XY;
for i = 2:1:size(width,1)
    idx = findCloset(libDist,width(i));
    if (width(i) >= (0.9 * libDist(idx))) && (width(i) <= (1.1 * libDist(idx)))
        model(i).XY = shape(idx).XY;
    else
        model(i).XY = model(i-1).XY;
    end
end
model(size(width,1)+1).XY = shape(end).XY;
% modelling
model11(model);
for d = 1:size(model,2)%% glitch
    tempA = smooth(model(d).XY);
    tempB = reshape(tempA,round(size(tempA,1)/2),2);
    smoothTest(d).XY = tempB;
    smoothTest(d).XY(1,1:2) = smoothTest(d).XY(end,1:2);
end
model11(smoothTest);
out.shape = shape;
out.model = model;
out.smooth = smoothTest;
out.width = width;
%save dataAlltoothpasteAD.mat;
end

function [idx]=findCloset(lib,val)
tmp = abs(lib-val);
[~, idx] = min(tmp); %index of closest value
%closest = lib(idx); %closest value
end
function []= model11(model)
reconX = [];reconY = [];reconZ = [];
shape = model;
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
fill3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','c')
plot3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','k-');
fill3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','c');
axis tight;
xlabel('Xaxis');ylabel('Yaxis');zlabel('height');
end
