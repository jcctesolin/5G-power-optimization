clc; clear; close all;

% Definição de N conjuntos
A{1} = {[1 0 0 0] [1 0 0 1] [1 1 0 0] [1 1 0 1]};
A{2} = {[0 0 1 0] [0 1 1 0] [0 0 1 1] [0 1 1 1]};
A{3} = {[0 1 0 0], [1 0 1 0], [1 1 0 0], [1 1 1 0]};


T = combinations(A{:});
V = T{:,:};
disp(T)

% T = combinations(A,B,C);
% V = T{:,:}; 
%disp(T)

%comb={};

% function z = comb(x,y)
%     for i=1:length(x)
%         for j=1:length(y)
%             z{i,j}= {x{i} y{j}};
%         end
%     end
% end
% 
% z = comb(A{1},A{2});
% disp(z);