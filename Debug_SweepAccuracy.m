%Sweep codes to analyze performance

%% First check bubble/dew point curves

%Grid of points to call
Ts = linspace(280, 450, 10); %K
NPs = 50;

%Compute results for each point
TPxys = zeros(size(Ts,1)*NPs, 4);

XNH3H2O('INIT', 'PropLib.dat');
ind = 1;
for j = 1:length(Ts)
    %First find P_min
    F_pmin = @(P) 0.0001 - XNH3H2O('TPX_Q', Ts(j), P, 0.0001);
    F_pmax = @(P) 0.9999 - XNH3H2O('TPX_Q', Ts(j), P, 0.9999);

    P_min = fzero(F_pmin, [10 5E7]);
    P_max = fzero(F_pmax, [10 5E7]);
    Ps = linspace(P_min, P_max, NPs);
    
    for k = 1:NPs        
        x = XNH3H2O('TPQ_X', Ts(j), Ps(k), 1E-6);
        y = XNH3H2O('TPQ_X', Ts(j), Ps(k), 1 - 1E-6);
        TPxys(ind,:) = [Ts(j) Ps(k) x y];
        ind = ind + 1;
    end
    
    %Quick check
    figure(1); clf; hold on;
    plot(TPxys(ind-NPs:ind-1, 3), TPxys(ind-NPs:ind-1, 2) );
    plot(TPxys(ind-NPs:ind-1, 4), TPxys(ind-NPs:ind-1, 2) );
    drawnow;
    pause(2);
    
end
XNH3H2O('END');


%% TPQ -> XY
Ts = linspace(270, 450, 10); %K
Ps = logspace(log10(1E3), log10(5E6), 10); %Pa
Qs = linspace(0.01, 0.99, 10);

[TTs, PPs, QQs] = meshgrid(Ts, Ps, Qs);
TTs = permute(TTs, [3 1 2]);
PPs = permute(PPs, [2 1 3]);
QQs = permute(QQs, [3 2 1]);

TPQs = [TTs(:) PPs(:) QQs(:)];
Xs = zeros(size(TPQs, 1), 1);
XNH3H2O('INIT', 'PropLib.dat');

for j = 1:length(Xs)
    Xs(j) = XNH3H2O('TPQ_X', TPQs(j,1), TPQs(j,2), TPQs(j,3) );
end
XNH3H2O('END');
