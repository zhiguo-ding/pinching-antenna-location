clear
%close all

ct = 5000;
snr=10000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [0:5 : 30];
sigmanoise = 10^(-12);
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =4*D;
Rtar = 2;
Rvec = [1:0.1:2];
  


for mi = 1 : length(Rvec) 
    Rtar = Rvec(mi);
    eps = sigmanoise/eta*(exp(M*Rtar)-1);
    %snr = 10^((snrvec(mi)-30)/10); 
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0; sum6=0;
    r_all = zeros(M,ct);
    r_all_pin = zeros(M,ct);
    r_all_pin2 = zeros(M,ct);
    for i = 1 : ct    
        loc = zeros(M,2);
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length        
        loc(:,2) = D*rand(M,1)-D/2; %width,    
         %with cluster
        loc(:,1) = D_leng/8*rand(M,1)-D_leng/4; %length 

        %conventional antennas OMA
        dall_conv = loc(:,1).^2+loc(:,2).^2+height^2;
        P_conv_temp(i) = sum(eps*dall_conv);

        %conventional antenna NOMA
        %pure NOMA
        d1 = dall_conv(1); d2 = dall_conv(2);
        P1conv_noma = sigmanoise*d1/eta*(exp(Rtar)-1);
        P2conv_noma = (exp(Rtar)-1)*((exp(Rtar)-1)*sigmanoise*d1+sigmanoise*d2)/eta;
        P_conv_noma_temp(i) = P1conv_noma+P2conv_noma;
        
        %pinching antennas OMA 
        xopt_sumpower = sum(loc(:,1))/M;
        dall_sumpower = [(xopt_sumpower-loc(:,1)).^2+loc(:,2).^2+height^2];
        P_pin_OMA_temp(i) = sum(eps*dall_sumpower);

        %exhaustive search of all possible positions
        stepx = D_leng/100;
        xvec = [-D_leng/2:stepx:D_leng/2];
        Pnoma_temp(i) = inf;
        Phnoma_temp(i) = inf;
        for ix = 1 : length(xvec)
            xpin = xvec(ix);
            
            %two distances
            dall = [(xpin-loc(:,1)).^2+loc(:,2).^2+height^2];

            d2 = max(dall);%weak user
            d1 = min(dall);%strong user

            %pure NOMA
            P1temp = sigmanoise*d1/eta*(exp(Rtar)-1);
            P2temp = (exp(Rtar)-1)*((exp(Rtar)-1)*sigmanoise*d1+sigmanoise*d2)/eta;
            if P1temp+P2temp<Pnoma_temp(i)
                Pnoma_temp(i)=P1temp+P2temp;
            end
 
            
        end

        % %use optimization solver        
        % A = []; % No other constraints
        % b = [];
        % Aeq = [];%
        % beq = [];%zeros(M+1,1); 
        % lb = [];
        % ub = [];
        % xsolver0 =  ones(3,1); %initilization 
        % options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        % epsnoma1 = sigmanoise/eta*(exp(Rtar)-1);
        % taum1 = epsnoma1*(loc(:,2).^2+height^2);
        % xsolver1 = fmincon(@(xsolver1) sum(xsolver1(1:2)),xsolver0,A,b,Aeq,beq,lb,ub,...
        %     @(xsolver1) myconsnoma(xsolver1,epsnoma1, eta, sigmanoise,loc,taum1 , D_leng)  ,options);
        % loc_temp = [loc(2,:);loc(1,:)];%swap SIC order
        % epsnoma2 = sigmanoise/eta*(exp(Rtar)-1);
        % taum2 = epsnoma2*(loc_temp(:,2).^2+height^2);
        % xsolver2 = fmincon(@(xsolver2) sum(xsolver2(1:2)),xsolver0,A,b,Aeq,beq,lb,ub,...
        %     @(xsolver2) myconsnoma(xsolver2,epsnoma2, eta, sigmanoise,loc_temp,taum2 , D_leng)  ,options);
        % 
        % Psolver_temp(i) = min(sum(xsolver1(1:2)), sum(xsolver2(1:2)));
        % xsolver = [xsolver1;xsolver2];
        % % one SIC order
        %     xpin = xsolver(1,3);            
        %     %two distances
        %     dall = [(xpin-loc(:,1)).^2+loc(:,2).^2+height^2]; 
        %     %rate
        %     P1temp = xsolver(1,1);P2temp = xsolver(1,2);
        %     if P1temp+P2temp<Pnoma_temp(i)
        %         Pnoma_temp(i)=P1temp+P2temp;
        %     end

        %analysis pure noma    
                
        %%%%%%channel order%%%%
        if abs(loc(1,2))>abs(loc(2,2)) %ensure user 1 is strong user
            a = loc(1,:);
            loc(1,:) = loc(2,:);
            loc(2,:) = a;
        end
        epsnoma = sigmanoise/eta*(exp(Rtar)-1);
        taum = epsnoma*(loc(:,2).^2+height^2);
        xop = loc(2,1)/(epsnoma*eta/sigmanoise+2) ...
            +loc(1,1)*(1+epsnoma*eta/sigmanoise)/(epsnoma*eta/sigmanoise+2);
        P1 = epsnoma*(xop-loc(1,1))^2+taum(1);
        P2 = epsnoma*eta/sigmanoise*P1+epsnoma*(xop-loc(2,1))^2+taum(2);
        Pnoma_temp_ana(i) = P1 + P2;

        %TEST
        %hnoma
 
    
     end
  
    Pconv_sim(mi) = sum(P_conv_temp)/ct;
    Pnoma_sim(mi) = sum(Pnoma_temp)/ct;
    Pnoma_ana(mi) = sum(Pnoma_temp_ana)/ct; 
    P_conv_noma(mi) = sum(P_conv_noma_temp)/ct;
    P_pin_OMA(mi) = sum(P_pin_OMA_temp)/ct;
    %Psolver(mi) = sum(Psolver_temp)/ct;
    
    
end
% plot(Rvec, Pconv_sim,Rvec, Pnoma_sim, Rvec, Psolver, Rvec,Pnoma_ana )
  plot(Rvec, Pconv_sim,Rvec, P_conv_noma,Rvec, P_pin_OMA,Rvec, Pnoma_sim, Rvec,Pnoma_ana )

 function [c,ceq] = myconsnoma(xsolver,epsnoma, eta, sigmanoise,loc,taum , D_leng)  

P1 = xsolver(1);
P2 = xsolver(2);
x = xsolver(3);
c(1,1) = -P2 + epsnoma*eta/sigmanoise*P1+epsnoma*(x-loc(2,1))^2+taum(2);
c(2,1) = -P2 + epsnoma*eta/sigmanoise*P1+epsnoma*(x-loc(1,1))^2+taum(1);
c(3,1) = -P1 + epsnoma*(x-loc(1,1))^2+taum(1);
c = [c;-P1;-P2];
c = [c;-D_leng/2-x;x-D_leng/2];
ceq = [];
end