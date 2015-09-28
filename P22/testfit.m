clear all;close all;
load testdatafit.mat;

% newSnap1 = rotateLine2(edgesnap(matching,1:2));
newSnap1 = edgesnap(matching,1:2);
x2 = newSnap1(:,1);
%x2 = x2(1:1:ceil(size(x2,1)/2));
y2 = newSnap1(:,2);
%y2 = y2(1:1:ceil(size(y2,1)/2));
% x = linspace(0,4*pi,10);
% y = sin(x);
% p = polyfit(x,y,7);
% 
% x1 = linspace(0,4*pi);
% y1 = polyval(p,x1);
% figure;
% plot(x,y,'o')
% hold on
% plot(x1,y1)
% hold off

%x2 = [2 4 6 7 9 11 14 22 ];
%x21 = [2 4 6];
%y2 = sin(x2);
name = 'poly';
ext = '.png';
a = 1;
while a <=10
p = polyfit(x2,y2,a);
x22 = x2;
%tInt = [-2:2:20];
y22 = polyval(p,x2);
fig = figure;
plot(x2,y2,'o');
hold on;
plot(x22,y22,'-r');
cat = strcat(name,num2str(a),ext);
saveas(fig,cat);
close all;
a = a+1;
end