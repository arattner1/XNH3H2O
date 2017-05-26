%For checking proplib validity
%Set the GetTPX function to always make a non-DB initial guess!

%Get values from the library
D = load('PropLib3.dat');
T_lib = D(:,1);
P_lib = D(:,2);
X_lib = D(:,3);
Q_lib = D(:,4);
H_lib = D(:,7);

%Make H error
Q_err = zeros(size(H_lib));
D_new = zeros(size(D));
%Loop through and calculate H error
XNH3H2O('INIT', 'PropLib2.dat');
for i = 1:length(T_lib)
    Q = XNH3H2O('TPX_Q', T_lib(i), P_lib(i), X_lib(i) );
    Q_err(i) = abs( (Q - Q_lib(i)) );
    %Make new D
    D_new(i, 1) = T_lib(i);
    D_new(i, 2) = P_lib(i);
    D_new(i, 3) = X_lib(i);
    D_new(i, 4) = Q;
    D_new(i, 5) = XNH3H2O('TPX_x', T_lib(i), P_lib(i), X_lib(i) );
    D_new(i, 6) = XNH3H2O('TPX_y', T_lib(i), P_lib(i), X_lib(i) );
    D_new(i, 7) = XNH3H2O('TPX_H', T_lib(i), P_lib(i), X_lib(i) );
    D_new(i, 8) = XNH3H2O('TPX_S', T_lib(i), P_lib(i), X_lib(i) );
    D_new(i, 9) = XNH3H2O('TPX_U', T_lib(i), P_lib(i), X_lib(i) );
    D_new(i, 10)= XNH3H2O('TPX_V', T_lib(i), P_lib(i), X_lib(i) );
end
XNH3H2O('END');