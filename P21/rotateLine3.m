function [newArea]=rotateLine3(img)
D1 = [img(:,2) img(:,1)];
[D1 ,~] = dpsimplify(D1,1);
% plot(D1(:,1),D1(:,2),'.r');
x_center = D1(1,2);
y_center = D1(1,1);
center = (repmat([x_center; y_center], 1, length(D1)));
theta = pi;
v = [D1(:,1)';D1(:,2)'];
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
last = repmat([max(D1(:,1))-min(D1(:,1));0],1,length(D1));
%vo = R*(v - center) + center + last;
vo = R*(v - center) + center +last ;
x_rotated = vo(1,:);
y_rotated = vo(2,:);

shiftX = y_center - x_rotated(1);
shiftY = x_center - y_rotated(end);
x_rotated = x_rotated + shiftX;
y_rotated = y_rotated + shiftY;

%hold on;plot(D1(:,2), D1(:,1), 'k-', x_rotated, y_rotated, 'r-', x_center, y_center, 'bo');

newArea = [D1 ; x_rotated' y_rotated'];
newArea = [newArea ;D1(1,1) D1(1,2)];
newArea = floor(newArea);
newArea = [newArea(:,2),newArea(:,1)];
%[newArea,~] = dpsimplify(newArea,1);
%D1 = newArea;
end