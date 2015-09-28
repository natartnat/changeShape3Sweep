function [width] = widthBasic(cropimage)
for a = size(cropimage,1)
    candidate(a,:) = find(cropimage(a,:) == searchFor);
    width(a,:) = candidate(a,end)-candidate(a,1);
    shape(a).XY;
end