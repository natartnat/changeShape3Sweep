function [pLayer] = B3swLoft(shape,width,height,cont)
%% get data to preparing resize profile
reconX = [];
reconY = [];
reconZ = [];

for i= 1:-1:((-1)*height)
    twoDresize = imresize(shape,[width width]);
    [yy , xx] = find(shape >= 1);
    %%
    recon = shape;
    reconX = [reconX ;recon(:,1)'];
    reconY = [reconY ;recon(:,2)'];
    reconZ = [reconZ ;(ones(1,size(reconX,2))*i)];    
end
    %pLayer(i).objLayer = [reconX,reconY,reconZ];
    pLayer = 1;
   hSurface = surf(reconX,reconY,reconZ,'EdgeColor','k','MeshStyle','column','FaceColor','interp','FaceLighting','gouraud');
set(hSurface,'FaceColor','c');
hold on;
plot3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','k-');
fill3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','c')
plot3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','k-');
fill3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','c');
end
