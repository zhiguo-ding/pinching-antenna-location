clear
%close all

ct = 50000;
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
Rtar = 3;
eps = sigmanoise/eta*(exp(M*Rtar)-1);  


parfor mi = 1: length(snrvec) 
    snr = 10^((snrvec(mi)-30)/10); 
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0; sum6=0; 
    Rate_conv_maxmin_temp = zeros(ct,1);
    Rate_conv_temp = zeros(ct,1);
    Rate_pin_sum_temp_sim = zeros(ct,1);
    Rate_pin_sum_temp_cub = zeros(ct,1);
    for i = 1 : ct    
        loc = zeros(M,2);
        %without cluster
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length  
        loc(:,2) = D*rand(M,1)-D/2; %width,

        % %with cluster
        % loc(:,2) = D/8*rand(M,1)-D/4; %width,
        loc(:,1) = D_leng/8*rand(M,1)-D_leng/4; %length 
        dall_conv = loc(:,1).^2+loc(:,2).^2+height^2; %tau_mx             

        % %conventional antennas - max-min
        % 
        % P_conv_temp = dall_conv./(sum(dall_conv))*snr;
        % Rate_conv_maxmin_temp(i) = sum(1/M*log(1+eta*P_conv_temp/sigmanoise./dall_conv));
        
        %conventional antenna - max sum rate
        %solution 1: user 2's constraint is met
        xtemp = eps*dall_conv(2);
        sol1 = [snr-xtemp xtemp];
        con1 = min([sol1+sol1/1000>=0  sol1+sol1/1000>=eps*dall_conv' snr+snr/1000>=sum(sol1)]);
        lam2 = 1/(snr-eps*dall_conv(2)+sigmanoise*dall_conv(1)/eta)...
            -1/(eps*dall_conv(2)+sigmanoise*dall_conv(2)/eta);
        ratez1 = sum(1/M*log(1+eta*sol1/sigmanoise./dall_conv'));
        %solution 2: user 1's constraint is met
        xtemp = eps*dall_conv(1);
        sol2 = [xtemp snr-xtemp];
        lam1 = 1/(snr-eps*dall_conv(1)+sigmanoise*dall_conv(2)/eta)...
            -1/(eps*dall_conv(1)+sigmanoise*dall_conv(1)/eta);
        con2 = min([sol2+sol2/1000>=0  sol2+sol2/1000>=eps*dall_conv' snr+snr/100>=sum(sol2)]);
        ratez2 = sum(1/M*log(1+eta*sol2/sigmanoise./dall_conv'));
        %solution 3: none of the equality constraints is met
        l0 = 2/(snr+sigmanoise*sum(dall_conv)/eta);
        sol3 = [1/l0 - sigmanoise*dall_conv(1)/eta 1/l0 - sigmanoise*dall_conv(2)/eta];
        con3 = min([sol3+sol3/1000>=0  sol3+sol3/1000>=eps*dall_conv' snr+snr/1000>=sum(sol3)]);
        ratez3 = sum(1/M*log(1+eta*sol3/sigmanoise./dall_conv'));
        if lam2>=0
            Rate_conv_temp(i) = ratez1;
        elseif lam1>=0
            Rate_conv_temp(i) = ratez2;
        else
            Rate_conv_temp(i) = ratez3;
        end
        %Rate_conv_temp(i) = max([ratez1*con1 ratez2*con2 ratez3*con3 ]);          

  

  
    end
    
  
    Rate_pin_sum_sim(mi) = sum(Rate_pin_sum_temp_sim)/ct; 
    Rate_pin_sum_cub(mi) = sum(Rate_pin_sum_temp_cub)/ct; 
    Rate_conv(mi) = sum(Rate_conv_temp)/ct;
    %Rate_conv_maxmin(mi) = sum(Rate_conv_maxmin_temp)/ct;
    
      
   
end
% plot(snrvec, Rate_conv_maxmin, snrvec, Rate_conv,snrvec, Rate_pin_sum_sim,snrvec, Rate_pin_sum_cub)
 plot(  snrvec, Rate_conv,snrvec, Rate_pin_sum_sim,snrvec, Rate_pin_sum_cub)
% pin = [loc [xopt;0] dall_opt/100 [index_opt;0]]