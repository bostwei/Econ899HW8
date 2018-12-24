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
k_space = [linspace(0.0001,K_ss,20),linspace(K_ss+0.01,14.99999,10)]'; 
K_space = [linspace(0.0001,K_ss * N,20),linspace(K_ss * N+0.01,14.99999 * N,10)]'; 

kk_space = linspace(0.1,14.99,2000)';
% construct a space with K and Z
[K_temp,e_emp,z_temp] = ndgrid(K_space,[1;0],[1;0]);
Kez_space = [K_temp(:),e_emp(:),z_temp(:)];

% calculate the future aggregate capital
KK = exp(a0 + a1 * log(Kez_space(:,1))) .* Kez_space(:,3) + exp(b0 + b1 * log(Kez_space(:,1))).* (1-Kez_space(:,3));

% calculate the accumulate labor choice
L = Kez_space(:,3)* L_g .* N +  (1 - Kez_space(:,3)) * L_b .* N;

% calclulate the technology shock
z = Kez_space(:,3) * zt(1) + (1 - Kez_space(:,3)) * zt(2); 

% calculate wage
w = (1 - alpha) .* z .*(Kez_space(:,1)./L).^ alpha;

 % calculate the interst rate
r = alpha .* z .*(Kez_space(:,1)./L).^ (1 - alpha); 

% the following will divide into 4 group
% -column 1 e = 1, z=1
% -column 2 e = 0, z=1
% -column 3 e = 1, z=0
% -column 4 e = 0, z=0
KK_temp = reshape(KK,size(k_space,1),2,2);
L_temp = reshape(L,size(k_space,1),2,2);
w_temp = reshape(w,size(k_space,1),2,2);
r_temp = reshape(r,size(k_space,1),2,2);
z_temp = reshape(z,size(k_space,1),2,2);

%% ----------------- generate person's choice given(k,z,e) ---------------
 inputDec.K_space =K_space;
 inputDec.kk_space = kk_space;

 inputDec.r_temp = r_temp;
 inputDec.K_temp = K_temp; 
% Note here the column in the dimention represent the following
% -column(:, 1,1) e = 1, z=1
% -column(:, 1,2) e = 0, z=1
% -column(:, 2,1) e = 1, z=0
% -column(:, 2,2) e = 0, z=0 
 inputDec.w_temp = w_temp;
 inputDec.e_emp = e_emp; % employment status
% Note here the column in the dimention represent the following
% -column(:, 1,1) e = 1, z=1
% -column(:, 1,2) e = 0, z=1
% -column(:, 2,1) e = 1, z=0
% -column(:, 2,2) e = 0, z=0 
 inputDec.e_bar = e_bar;
 inputDec.delta = delta;
 inputDec.beta = beta;
 inputDec.N = N;
 
 tic
 outDec = Dec(inputDec);
 toc



%%

%------------ calcuate the personal choice of capital-----------------
% initial the asset choice for different people and 

ind_k = zeros(N,T); % individual asset choice.

% extend the Z to each individual
Z_ext = ones(N,1) * Z;

% start with the initial asset holding to be the K_ss
ind_k(:,1) = K_ss * ones(N,1);

% initiate totol asset holding
K_sum = zeros(1,T);
K_sum(1) = K_ss;

% calculate the aggregate labor supply
L_sum = (L_g .* Z + L_b .*(1- Z)).*N;

% initiate wage sequence
w_seq = zeros(1,T);
w_seq(1) = w;

% initiate rate sequence
r_seq = zeros(1,T);
r_seq(1) = r;

%%
%Next period's asset holding depend by the interpolation of (k,e,z)
for t = 2: 2
      % calculate wage
    w_seq(t) = (1 - alpha) .* Z_ext(:,t-1) .*(K_sum(:,t-1)./L_sum(:,t)).^ alpha;

     % calculate the interst rate
    r_seq(t) = alpha .*  Z_ext(:,t-1) .*(K_sum(:,t-1)./L_sum(:,t)).^ (1 - alpha); 
 
 
    
    % determined the asset holding next period
    % e = 1, z = 1
        k_1 = interp1(input.k(:,1,1),dec_k(:,1,1),...
           ind_k(:,t-1),'spline').* employed(:,t-1) .*Z_ext(:,t-1);
    % e = 0, z = 1 
        k_1 = k_1 + interp1(input.k(:,1,2),dec_k(:,1,2),...
       ind_k(:,t-1),'spline').*(1 - employed(:,t-1)) .*Z_ext(:,t-1);
    % e = 1, z = 0 
        k_1 = k_1 + interp1(input.k(:,2,1),dec_k(:,2,1),...
       ind_k(:,t-1),'spline').*employed(:,t-1) .*(1 - Z_ext(:,t-1));
    % e = 0, z = 0 
        k_1 = k_1 + interp1(input.k(:,2,2),dec_k(:,2,2),...
        ind_k(:,t-1),'spline').*(1-employed(:,t-1)) .*(1 - Z_ext(:,t-1));
    
    ind_k(:,t) = k_1;
    
    % sum of aggregate asst 
    K_sum(:,t) = sum(ind_k(:,t),1);
    
        
    

    

end

