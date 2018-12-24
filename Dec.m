function [outputArg1] = Dec(inputArg)
%DEC calcualte the decision matrix given k,e,K,z
%   Detailed explanation goes here

K_space = inputArg.K_space;
kk_space = inputArg.kk_space;
N = inputArg.N;

v0 = ones(size(K_space,1),2,2)*100; % the v0(k,e,z)

uinput.r = inputArg.r_temp;
uinput.k = inputArg.K_temp./N; 
% Note here the column in the dimention represent the following
% -column(:, 1,1) e = 1, z=1
% -column(:, 1,2) e = 0, z=1
% -column(:, 2,1) e = 1, z=0
% -column(:, 2,2) e = 0, z=0 
uinput.w = inputArg.w_temp;
uinput.e = inputArg.e_emp; % employment status
% Note here the column in the dimention represent the following
% -column(:, 1,1) e = 1, z=1
% -column(:, 1,2) e = 0, z=1
% -column(:, 2,1) e = 1, z=0
% -column(:, 2,2) e = 0, z=0 
uinput.e_bar = inputArg.e_bar;
uinput.delta = inputArg.delta;

uinput.beta = inputArg.beta ; 

tol = 10;
Iter = 0;
while tol > 1
    
uinput.v0 = v0;
% calculate the desition choice given v0 and kk_space
out = u(uinput,kk_space);

% update v1 and v0;
v1 = out.v1;
dec = out.dec_k;
dec_k = kk_space(dec);

tol = max(max(max(abs(v1-v0))));

v0 = v1;
Iter = Iter+1;
fprintf('This is %d iteration, the distance is %.4f .\n',Iter,tol);
end

outputArg1.v1 = v1;
outputArg1.v0 = v0;
outputArg1.dec = dec;
outputArg1.dec_k = dec_k;
end

