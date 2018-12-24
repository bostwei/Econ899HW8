%-------------------------------------------------------------------------
% Econ 899 Hw8
% Shiyan Wei
% 12/22/2018
%-------------------------------------------------------------------------

%-----------------------------------------------------------------------
% This script is to repicate the Aiyagari (1994) paper with aggregate
% uncertainty setting using Krusell and Smith (1998) techiques.
%-----------------------------------------------------------------------

clear
clc

%------------ parameter initiated ------------------------------

alpha = 0.36; % preference substitude 
beta = 0.99; % discount rate
delta = 0.025; % capital depreciation rate
epsilon= [0;1]; % the employment opportunity state
zt = [1.01;0.99]; % the technology space
e_bar = 0.3271; % labor efficiency per unity of work time

N = 5000; % number of simulation people 
T = 11000; % period of time simulated 

%-------------------- generate the transition matrix -----------------
 run('transmat.m')
 
 
 % generate conditional probabilities
% - row is the current state
% - coloumn is the future state
 pr_con = [pgg11, pgb11, pgg10, pgb10;
       pbg11, pbb11, pbg10, pbb10;
       pgg01, pgb01, pgg00, pgb00;
       pbg01, pbb01, pbg00, pbb00];
   
% %  pi_gg = [pgg11, pgg01;
% %           pgg10, pgg00];
%  pi_gg = [ pgg11,  1-pgg00;
%           1-pgg11, pgg00];
% %  pi_gb = [pgb11, pgb01;
% %           pgb10, pgb00];
%  pi_gb = [ pgb11,  1-pgb00;
%           1-pgb11, pgb00];

% after transition
pr_con = pr_con';
% - column is the current state
% - row is the future state

% labor supply given status
L_g = 0.96;
L_b = 0.9;


% Initial K
K_ss = 5.7163;


%----------------- determined the state status ----------------------
Z = zeros(1,T); % the vectore indicate the good or bad state shocks 
% z = 1; good shock
% z = 0; bad shock
random = rand(1,T); % generate the randome shock  

% detect whether current state will change the state
change_State = double(random >0.875); % =1 change the current state to be different than the previous state
Z(1) = 1;

%
for i = 2: T
    if change_State(i) == 1 % if change state need to change
        if Z(i-1) == 1
           Z(i) = 0;
        elseif Z(i-1) ==0
           Z(i) = 1; 
        end
    else % if change state does not need to change
        Z(i) = Z(i-1);
    end
end

% ------------------- determined the emplotement status------------------ 
 
employed = zeros(N,T); % =1 if employed

% use random number to dertermined the employment status
random_emp = rand(N,T);

% update the first period unemployment status
% the first period is good status, therefore the aggregate unemployment
% rate is 4%, then for those random_emp > 0.96 are umemployed.

employed(:,1) = double(random_emp(:,1) <= 0.96);

% update the unemployment status for t>1
tic
for t = 2: T
    for n = 1:N
        % if the prevous state is good state and current state is also good state
        if Z(t-1) ==1 && Z(t) ==1 
            if employed(n,t-1) == 1 % if previous state is employed
               employed(n,t) = double(random_emp(n,t) <= pgg11); % then the current period of employed probability is rand < pgg11
            else % if the previous state is unemployed
               employed(n,t) = double(random_emp(n,t) <= pgg10); % then the current period of employed the future employment probability is rand <= pgg10. Note that the notation of pgg10 accually should write as pgg01  
            end
            
        % if the prevous state is good state and current state is bad    
        elseif Z(t-1) ==1 && Z(t) ==0 
            if employed(n,t-1) == 1 % if previous state is employed
               employed(n,t) = double(random_emp(n,t) <= pgb11); % then the current period of employed probability is rand < pgb11
            else % if the previous state is unemployed
               employed(n,t) = double(random_emp(n,t) <= pgb10); % then the current period of employed the future employment probability is rand <= pgb10. Note that the notation of pgg10 accually should write as pgg01  
            end   
            
        % if the prevous state is bad state and current state is good    
        elseif Z(t-1) ==0 && Z(t) ==1 
            if employed(n,t-1) == 1 % if previous state is employed
               employed(n,t) = double(random_emp(n,t) <= pbg11); % then the current period of employed probability is rand < pbg11
            else % if the previous state is unemployed
               employed(n,t) = double(random_emp(n,t) <= pbg10); % then the current period of employed the future employment probability is rand <= pbg10. Note that the notation of pgg10 accually should write as pgg01  
            end   
        % if the prevous state is bad state and current state is good 
        else 
            if employed(n,t-1) == 1 % if previous state is employed
               employed(n,t) = double(random_emp(n,t) <= pbg11); % then the current period of employed probability is rand < pbg11
            else % if the previous state is unemployed
               employed(n,t) = double(random_emp(n,t) <= pbg10); % then the current period of employed the future employment probability is rand <= pbg10. Note that the notation of pgg10 accually should write as pgg01  
            end
            
        end
            
    end
end
toc




%% Step 3 conjuncture the policy function

%initial guess the paramter in log form
a0 = 0.095;
b0 = 0.085;
a1 = 0.999;
b1 = 0.999;

% interpolate 10 point in k = K = [0,15)
k_space = linspace(0.1,14.99,20)'; 
K_space = linspace(0.1,14.99,20)';

kk_space = linspace(0.1,14.99,2000)';
% construct a space with K and Z
[K_temp,e_emp,z_temp] = ndgrid(K_space,[1;0],[1;0]);
Kez_space = [K_temp(:),e_emp(:),z_temp(:)];

% calculate the future aggregate capital
KK = exp(a0 + a1 * log(Kez_space(:,1))) .* Kez_space(:,2) + exp(b0 + b1 * log(Kez_space(:,1))).* (1-Kez_space(:,2));

% calculate the accumulate labor choice
L = Kez_space(:,2)* L_g  +  (1 - Kez_space(:,2)) * L_b;

% calclulate the technology shock
z = Kez_space(:,2) * zt(1) + (1 - Kez_space(:,2)) * zt(2); 

% calculate wage
w = (1 - alpha) .* z .*(Kez_space(:,1)./L).^ alpha;

 % calculate the interst rate
r = alpha .* z .*(Kez_space(:,1)./L).^ (1 - alpha); 


KK = reshape(KK,2,20)';
L = reshape(L,2,20)';
w = reshape(w,2,20)';
r = reshape(r,2,20)';
z = reshape(z,2,20)';

%----------------- generate person's choice given(z,e) ---------------

% calculate the utility given all these parameter
v0 = zeros(size(K_space,1),2,2); % the v0(k,e,z)


%% assembly all the input requried to calculate utility
input.r = r;
input.k = K_temp';
input.w = w;
input.e = z_temp'; % employment status

input.e_bar = e_bar;
input.delta = delta;

input.beta = beta ; 

tol = 10;
Iter = 0;
while tol > 0.001
    
input.v0 = v0;
% calculate the desition choice given v0 and kk_space
out = u(input,kk_space);

% update v1 and v0;
v1 = out.v1;
dec = out.dec_k;
dec_k = kk_space(dec);

tol = max(max(abs(v1-v0)));

v0 = v1;
Iter = Iter+1;
fprintf('This is %d iteration, the distance is %.4f .\n',Iter,tol);
end




%%

%------------ calcuate the accumulative capital and Labor choice-----------------

% % calculate the accumulative capital 
% K = zeros(1,T);
% K(1) = K_ss;
% 
% 
% 
% 
% L = zeros(1,T);
% z = zeros(1,T);
% w = zeros(1,T);
% r = zeros(1,T);
% 
% V = zeros(N,T); % preson n's capital choice given the state z
% for t = 1:T-1
%     % if current period is good 
%     if Z(t) == 1 
%        K(t+1) = exp(a0 + a1 * log(K(t)));
%     
%     % if the current period is bad
%     else
%        K(t+1) = exp(b0 + b1 * log(K(t)));
%     end
% 
% 
% % calculate the accumulate labor choice
% L(t) = Z(t)* L_g  +  (1 - Z(t)) * L_b;
% 
% % calclulate the technology shock
% z(t) = Z(t) * zt(1) + (1 - Z(t)) * zt(2); 
% 
% % calculate wage
% w(t)= (1 - alpha) .* z .*(K(t)./L(t)).^ alpha;
%  
%  % calculate the interst rate
% r(t) = alpha .* z .*(K(t)./L(t)).^ (1 - alpha); 
% 
% end