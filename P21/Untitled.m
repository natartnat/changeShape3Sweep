%% unused code 
% get width in data line
%     [mMin, imin] = min(xy1(:,1));
%     start_x = mMin;
%     start_y = xy1(imin,2);
%     
%     [mMax, imax] = max(xy1(:,1));
%     stop_x = mMax;
%     stop_y = xy1(imax,2);
%     
%     width(1) = mMax - mMin;

% for a = 2:divide
%     %width(a) = max(xy1,1) - min(xy1,1);
%     %search around + - 10%
%     oneTh = ceil(width(a-1)/10);
%     %temp = find(ff(ceil(start_y+1),:) <= start_x + oneTh & ff(ceil(start_y+1),:)>= start_x - oneTh);
%     temp = ff(ceil(start_x-oneTh):ceil(start_x-oneTh),ceil(start_y+1));
%     
%     if ~isempty(temp)
%         start_x = ff(start_y+1,temp(1)); %select first in 10%
%         start_y = start_y+1;
%     else
%         start_y = start_y+1;
%     end
%     
%     temp = find(ff(ceil(stop_y+1),:) <= stop_x + oneTh & ff(ceil(stop_y+1),:)>= stop_x - oneTh);
%     if ~isempty(temp)
%         stop_x = ff(stop_y+1,temp(1)); %select first in 10%
%         stop_y = stop_y+1;
%     else
%         stop_y = stop_y+1;
%     end
%     width(a) = stop_x - start_x;
% end
% width = width';
%model=B3swLoft(shape,width,divide);
%%
% |MONOSPACED TEXT|
%Specify initial contour interactively.
%mask = roipoly;
%figure, imshow(A);
%figure, imshow(mask)
% title('Initial MASK');

%%%Segment the image, specifying 200 iterations.
% maxIterations = 200; 
% bw = activecontour(I, mask, maxIterations, 'Chan-Vese');
%   
% % Display segmented image
% figure, imshow(bw)
% title('Segmented Image');