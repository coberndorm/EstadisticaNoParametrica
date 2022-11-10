clear
filename = 'return.txt';
[X,~] = importdata(filename);

%% Punto 1
X900 = X(end-899:end,:);

%