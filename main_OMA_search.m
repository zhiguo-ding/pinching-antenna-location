clear
%close all

ct = 50;
snr=10000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [0:5 : 20];
sigmanoise = 10^(-12);
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =3*D;
Rtar = 2;
Rvec = [2:1:6];
  


for mi = 1 : length(Rvec) 
    Rtar = Rvec(mi);
    eps = sigmanoise/eta*(exp(M*Rtar)-1);
    snr = 10^((snrvec(mi)-30)/10); 
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0; sum6=0;
    r_all = zeros(M,ct);
    r_all_pin = zeros(M,ct);
    r_all_pin2 = zeros(M,ct);
    for i = 1 : ct    
        loc = zeros(M,2);
        %without cluster
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length  
        loc(:,2) = D*rand(M,1)-D/2; %width,
        %with cluster
        %loc(:,2) = D/8*rand(M,1)-D/4; %width,
        loc(:,1) = D_leng/8*rand(M,1)-D_leng/4; %length 
             

        %conventional antennas
        dall_conv = loc(:,1).^2+loc(:,2).^2+height^2
        P_conv_temp(i) = sum(eps*dall_conv);

        %exhaustive search of all possible positions
        stepx = D_leng/100;
        xvec = [-D_leng/2:stepx:D_leng/2];
        P_sumpower_temp(i) = inf;
        Phnoma_temp(i) = inf;
        for ix = 1 : length(xvec)
            xpin = xvec(ix);
            
            %all distances
            dall = [(xpin-loc(:,1)).^2+loc(:,2).^2+height^2]; 
            %minimize sum power
            P_sumpower_tempx = sum(eps*dall);
            if P_sumpower_tempx<P_sumpower_temp(i)
                P_sumpower_temp(i)=P_sumpower_tempx;
            end 
        end 



        %analysis sum power 
        xopt_sumpower = sum(loc(:,1))/2;
        dall_sumpower = [(xopt_sumpower-loc(:,1)).^2+loc(:,2).^2+height^2];
        P_sumpower_anax(i) = sum(eps*dall_sumpower);
    end
    
  
    Pconv_sim(mi) = sum(P_conv_temp)/ct; 
    P_sumpower_sim(mi) = sum(P_sumpower_temp)/ct;
    P_sumpower_ana(mi) = sum(P_sumpower_anax)/ct;
     
   
end
 plot(Rvec, Pconv_sim,Rvec, P_sumpower_sim, Rvec,P_sumpower_ana)