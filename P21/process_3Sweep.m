function data = process_3Sweep(inx1,iny1,inx2,iny2,countD,I)
disMax=0;
objNum = 1;
data = [];
for count = 1:1:countD-1
%% ============= make flood fill 2D Base =================
temp2D = selectOBJ(I,[round(inx1(count));round(iny1(count))]);
%     s  = regionprops(temp2D,'ConvexImage');
%     imshow(s.ConvexImage);
    
    %for k = 1:1:size(s.ConvexImage,1) 
        %hold on,plot(s.ConvexImage(k,1),s.ConvexImage(k,2),'b.') % TR RT
    %end
%figure,imshow(temp2D);
%% ================ 2D-BASAE construction =================
%[ twoD ] = twoDBaseProfile(temp2D);
[ twoD ] = twoDBaseProfile(temp2D);
%figure,imshow(twoD);

tempBody = selectOBJ(I,[round(inx2(count));round(iny2(count))]);
%figure,imshow(tempBody);
%% ================ MOD for 2 side Object =================
tempBody = select2Side( tempBody , temp2D , I);
%figure,imshow(twoD);
%figure,imshow(tempBody);
[height, profWidth] = getProperty(tempBody); %% length to place 2D
%[height, profWidth] = getPropertyV4(tempBody,temp2D); %% length to place 2D
%[height, profWidth] = getpropertyV2(tempBody); %% length to place 2D

if (iny2-iny1) > 0
    direction = -1;
    %['DOWN'];
else
    direction = 1;
    %['UP'];
end
% Basic3sw( twoD,distance,direction);
data = Basic3sw( twoD,profWidth,height,disMax,direction,objNum,data);
%
disMax = disMax + height;
objNum = objNum + 1;
end
axis tight;

%['END'];