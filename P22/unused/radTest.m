clear all;close all;
load testdatafit.mat;
newSnap1 = edgesnap(matching,1:2);
X = newSnap1(:,1);
T = newSnap1(:,2);
%X = -1:.1:1;
%T = [-.9602 -.5770 -.0729  .3771  .6405  .6600  .4609 ...
%     .1336 -.2013 -.4344 -.5000 -.3930 -.1647  .0988 ...
%     .3072  .3960  .3449  .1816 -.0312 -.2189 -.3201];
plot(X,T,'+');
title('Training Vectors');
xlabel('Input Vector P');
ylabel('Target Vector T');

x = -3:.1:3;
a = radbas(x);
plot(x,a)
title('Radial Basis Transfer Function');
xlabel('Input p');
ylabel('Output a');

a2 = radbas(x-1.5);
a3 = radbas(x+2);
a4 = a + a2*1 + a3*0.5;
plot(x,a,'b-',x,a2,'b--',x,a3,'b--',x,a4,'m-')
title('Weighted Sum of Radial Basis Transfer Functions');
xlabel('Input p');
ylabel('Output a');

eg = 0.02; % sum-squared error goal
sc = 1;    % spread constant
net = newrb(X,T,eg,sc);

plot(X,T,'+');
xlabel('Input');

X = -1:.01:1;
Y = net(X);

hold on;
plot(X,Y,'-r');
hold off;
legend({'Target','Output'})