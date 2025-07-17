clear
%close all

ct = 500000;
snr=10000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [0:2 : 20];
sigmanoise = 10^(-12);
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =4*D;
Rtar = 2.5;
  
parfor mi = 1 : length(snrvec) 
    snr = 10^((snrvec(mi)-30)/10); 
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0; sum6=0;
    r_all = zeros(M,ct);
    r_all_pin = zeros(M,ct);
    r_all_pin2 = zeros(M,ct);
    for i = 1 : ct    
        loc = zeros(M,2);
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length        
        loc(:,2) = D*rand(M,1)-D/2; %width,    
        loc(:,1) = D_leng/8*rand(M,1)-D_leng/4; %length 

        %conventional antenna
        dist = sqrt(loc(:,1).^2+ loc(:,2).^2 +height^2);
        r_all(:,i) = log(1+eta.*snr./dist.^2/sigmanoise)/M;      
        if r_all(1,i)< Rtar
            sum1 = sum1 + 1;
        end 
        
        %pinching antenna
        %case III
        %wd_loss = 10.^(0.08*(abs(loc(:,1)+D_leng/2))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %case II
        %wd_loss = 10.^(0.08*(abs(loc(:,1)))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %wd_loss = 1; %no loss for now
        % individual
 
        
        %fixed pinching antennas
        xpin = sum(loc(:,1))/M;
        dist3 = sqrt((loc(:,1)-xpin).^2+ loc(:,2).^2 +height^2);
        r_all_pin2(:,i) = log(1+eta*snr./dist3.^2/sigmanoise )/M;   
        if r_all_pin2(1,i)< Rtar
            sum3 = sum3 + 1;
        end
 
    end
 
    %outage probability
    po_conv(mi) = sum1/ct; 
    po_pin(mi) = sum3/ct; 

    %analyticla results
    eps = sigmanoise/eta*(exp(M*Rtar)-1);
    if snr/eps-height^2<0
        pana(mi) = 1;
    else        
        theta1 = 4/D_leng^2*(snr/eps-height^2);
        theta4 =              D_leng^2/4*theta1;
        theta2 = min(D/2,sqrt(theta4)); 
        theta3 = min(theta2, sqrt(max(0,D_leng^2/4*(theta1-1))));
        
        y1=theta2;
        temp1 = y1+theta1*y1-4/3/D_leng^2*y1^3-4/D_leng*...
            (y1/2*sqrt(max(0,theta4-y1^2)) + theta4/2*asin(y1/sqrt(theta4)));

        y1=theta3;
        temp2 = y1+theta1*y1-4/3/D_leng^2*y1^3-4/D_leng*...
            (y1/2*sqrt(theta4-y1^2) + theta4/2*asin(y1/sqrt(theta4)));
        if theta1<0
            dfd=0;
        end
        pana(mi) = 2*(D/2-min(D/2,sqrt(theta1*D_leng^2/4)) )/D...
            +2/D*(temp1-temp2);   

    end
 
    
   
end
 
%semilogy(snrvec, po_conv, snrvec, po_1pin, snrvec, po_pin, snrvec,pana_old, snrvec,pana)
%semilogy(snrvec, po_conv, snrvec, po_pin, snrvec,pana)


plot(snrvec, (1-po_conv)*Rtar, snrvec, (1-po_pin)*Rtar, snrvec,(1-pana)*Rtar)
 