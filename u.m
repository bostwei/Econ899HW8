function output = u(input,kk)
%U function u calculate the utility function given k e  z; K,
%   Detailed explanation goes here
% Note here the column in the dimention represent the following
% -column(:, 1,1) e = 1, z=1
% -column(:, 1,2) e = 0, z=1
% -column(:, 2,1) e = 1, z=0
% -column(:, 2,2) e = 0, z=0 

r = input.r;
k = input.k;
w = input.w;
e = input.e; % employment status

% % reshape the variable to the form (k,e,z)
% r = reshape(r,size(k,1),2,2);
% w = reshape(w,size(k,1),2,2);


e_bar = input.e_bar;
delta = input.delta;
beta = input.beta; 



% value function before the interpolation
v0 = input.v0; % the v0(k,e,z)



% the is the form of Income(k,e,z)
% Income from employed good state
Income = r .* k + w .* e .* e_bar + (1 - delta) .* k; 

% the c is the form c(k,k',e,z)
c(:,:,1,1) = Income(:,1,1) - kk'; % consumption if e = 1, z = 1 
c(:,:,1,2) = Income(:,1,2) - kk'; % consumption if e = 1, z = 2 
c(:,:,2,1) = Income(:,2,1) - kk'; % consumption if e = 1, z = 1 
c(:,:,2,2) = Income(:,2,2) - kk'; % consumption if e = 1, z = 2 

c(c<0) = 0;

u = log(c);


% interpolation of the value function v0
kk_I = kk; % interpolation point of kk
v0_I(:,1,1) = interp1(k(:,1,1),v0(:,1,1),kk_I);
v0_I(:,1,2) = interp1(k(:,1,2),v0(:,1,1),kk_I);
v0_I(:,2,1) = interp1(k(:,2,1),v0(:,1,1),kk_I);
v0_I(:,2,2) = interp1(k(:,2,2),v0(:,1,1),kk_I);


% % calculate the wealth
% wealth for e = 1, z = 1
 wealth(:,:,1,1) = u(:,:,1,1) + beta .*  v0_I(:,1,1)';
 % wealth for e = 1, z = 1
 wealth(:,:,1,2) = u(:,:,1,2) + beta .*  v0_I(:,1,2)';
 % wealth for e = 1, z = 1
 wealth(:,:,2,1) = u(:,:,2,1) + beta .*  v0_I(:,2,1)';
 % wealth for e = 1, z = 1
 wealth(:,:,2,2) = u(:,:,2,2) + beta .*  v0_I(:,2,2)';


 
 % find the optimal choice of asset and value function
 [v1,dec_k] = max(wealth,[],2);
v1_temp = reshape(v1,size(k,1),2,2);
dec_k_temp = reshape(dec_k,size(k,1),2,2);
output.v1 = v1_temp;
output.dec_k = dec_k_temp;

end

