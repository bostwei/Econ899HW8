
%parameters of transition matrix:

durug=1.5;% duration of unemployment spell in good state
durub=2.5; % duration of unemployment spell in bad state

durgd=8.0; % duration of good state
durbd=8.0; % duration of bad state

unempg=0.04; % good time unemployment rate
unempb=0.1; % bad time unemployment rate



%transition probabilities of stata status
pgg = (durgd-1)/durgd;
pgb = 1 - (durbd-1)/durbd;
pbg = 1 - (durgd-1)/durgd;
pbb = (durbd-1)/durbd;


% transition of unemployment opportunity condition of status
pgg00 = (durug-1)/durug;
pbb00 = (durub-1)/durub;
pbg00 = 1.25*pbb00;
pgb00 = 0.75*pgg00;

pgg01 = (unempg - unempg*pgg00)/(1-unempg);
pbb01 = (unempb - unempb*pbb00)/(1-unempb);
pbg01 = (unempb - unempg*pbg00)/(1-unempg);
pgb01 = (unempg - unempb*pgb00)/(1-unempb);


pgg10 = 1 - (durug-1)/durug;
pbb10 = 1 - (durub-1)/durub;
pbg10 = 1 - 1.25*pbb00;
pgb10 = 1 - 0.75*pgg00;

pgg11 = 1 - (unempg - unempg*pgg00)/(1-unempg);
pbb11 = 1 - (unempb - unempb*pbb00)/(1-unempb);
pbg11 = 1 - (unempb - unempg*pbg00)/(1-unempg);
pgb11 = 1 - (unempg - unempb*pgb00)/(1-unempb);


%matrix
pr(1,1) = pgg*pgg11;
pr(2,1) = pbg*pbg11;
pr(3,1) = pgg*pgg01;
pr(4,1) = pbg*pbg01;
pr(1,2) = pgb*pgb11;
pr(2,2) = pbb*pbb11;
pr(3,2) = pgb*pgb01;
pr(4,2) = pbb*pbb01;
pr(1,3) = pgg*pgg10;
pr(2,3) = pbg*pbg10;
pr(3,3) = pgg*pgg00;
pr(4,3) = pbg*pbg00;
pr(1,4) = pgb*pgb10;
pr(2,4) = pbb*pbb10;
pr(3,4) = pgb*pgb00;
pr(4,4) = pbb*pbb00;



% pr = [pgg*pgg11, pgb*pgb11, pgg*pgg10, pgb*pgb10;
%      pbg*pbg11, pbb*pbb11, pbg*pbg10, pbb*pbb10;
%      pgg*pgg01, pgb*pgb01, pgg*pgg00, pgb*pgb00;
%      pbg*pbg01, pbb*pbb01, pbg*pbg00, pbb*pbb00];
% pr is the unconditional probabilities tranistion probabilities of unemployment state given different status(z,e)
% This is a transition Pi_zz'ee';




pr_tmp = pr^10000; % calculate the stationary probabilities transition


% generate conditional probabilities
 pr_con = [pgg11, pgb11, pgg10, pgb10;
       pbg11, pbb11, pbg10, pbb10;
       pgg01, pgb01, pgg00, pgb00;
       pbg01, pbb01, pbg00, pbb00];


for i=1:4
    pr_star(i) = pr_tmp(i,i); % stationary probabilities distribution of four state
end

% disp(['Stationary Distribution = [',num2str(pr_star),']'])
