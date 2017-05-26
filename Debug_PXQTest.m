%Debugging code to track down PXQ behavior error
%Alex Rattner, 2017-05-25
%Based on user feedback

%Restart code:
XNH3H2O('END');
XNH3H2O('INIT', 'PropLib2.dat');


%Inputs:
P = 850E3;
Xs = 0.90:0.001:0.999;
%Xs = [0.981];
Q  = 0.99;

%Out:
Hs = zeros(size(Xs));
Ts = zeros(size(Xs));

for j = 1:length(Hs)
    Hs(j) = XNH3H2O('PXQ_H', P, Xs(j), Q);
    
%    pause;
    
    Ts(j) = XNH3H2O('PXQ_T', P, Xs(j), Q);
    
%    pause;
end

figure(1); clf;
subplot(1,2,1);
plot(Xs, Hs);
subplot(1,2,2);
plot(Xs, Ts);