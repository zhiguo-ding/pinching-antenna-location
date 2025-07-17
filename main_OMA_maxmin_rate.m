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
  


parfor mi = 1 : length(snrvec) 
    snr = 10^((snrvec(mi)-30)/10); 
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0; sum6=0; 
    Rate_conv_temp = zeros(ct,1);
    Rate_pin_min_temp_sim = zeros(ct,1);
    Rate_conv_ana_temp = zeros(ct,1);
    for i = 1 : ct    
        loc = zeros(M,2);
        %without cluster
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length  
        loc(:,2) = D*rand(M,1)-D/2; %width,
        % %with cluster
        % %loc(:,2) = D/8*rand(M,1)-D/4; %width,
         loc(:,1) = D_leng/8*rand(M,1)-D_leng/4; %length 
             

        %conventional antennas
        dall_conv = loc(:,1).^2+loc(:,2).^2+height^2; %tau_mx
        P_conv_temp = dall_conv./(sum(dall_conv))*snr;
        Rate_conv_temp(i) = min(1/M*log(1+eta*P_conv_temp/sigmanoise./dall_conv));

        %exhaustive search of all possible positions
        stepx = D_leng/100;
        xvec = [-D_leng/2:stepx:D_leng/2];
        Rate_pin_min_temp_sim(i) = 0;
        for ix = 1 : length(xvec)
            xpin = xvec(ix);
            
            %all distances
            dall = [(xpin-loc(:,1)).^2+loc(:,2).^2+height^2]; 
            P_pin_temp = dall./(sum(dall))*snr;
            Rate_pin_temp = min(1/M*log(1+eta*P_pin_temp/sigmanoise./dall))
            if Rate_pin_temp>Rate_pin_min_temp_sim(i)
                Rate_pin_min_temp_sim(i)=Rate_pin_temp;
                xopt_search = xpin;
            end 
        end 



        %analysis sum power 
        xopt = sum(loc(:,1))/M;
        dall_pin_opt = (xopt-loc(:,1)).^2+loc(:,2).^2+height^2; %tau_mx
        P_pin_ana_temp = dall_pin_opt./(sum(dall_pin_opt))*snr;
        Rate_conv_ana_temp(i) = min(1/M*log(1+eta*P_pin_ana_temp/sigmanoise./dall_pin_opt));
    end
    
  
    Rate_pin_min_sim(mi) = sum(Rate_pin_min_temp_sim)/ct; 
    Rate_pin_min_ana(mi) = sum(Rate_conv_ana_temp)/ct; 
    Rate_conv(mi) = sum(Rate_conv_temp)/ct;
      
   
end
 plot(snrvec, Rate_conv,snrvec, Rate_pin_min_sim, snrvec,Rate_pin_min_ana)