         

 
        %exhaustive search of all possible positions
        stepx = D_leng/100;
        xvec = [-D_leng/2:stepx:D_leng/2];
        Rate_pin_sum_temp_sim(i) = 0;
        for ix = 1 : length(xvec)
            xpin = xvec(ix);
            
            %all distances
            dall = [(xpin-loc(:,1)).^2+loc(:,2).^2+height^2];   
            % objf(ix) = log(1/2+eta*(snr/2+sigmanoise*dall(2)/2/eta)/sigmanoise/dall(1))...
            %     +log(1/2+eta*(snr/2+sigmanoise*dall(1)/2/eta)/sigmanoise/dall(2));
            % objf(ix) = dall(1)/4/dall(2) + dall(2)/4/dall(1) + ...
            %     eta*snr/2/sigmanoise/dall(1) + eta*snr/2/sigmanoise/dall(2)...
            %     +eta^2*snr^2/4/sigmanoise^2/dall(1)/dall(2);           
            %objf(ix) = eta^2*snr^2/4/sigmanoise^2/dall(1)/dall(2);
            objf(ix) = dall(1)*dall(2);
        end  
        plot(xvec,objf)
        [xxxz,xxxxy] = min(objf);
        [xopt xvec(xxxxy)]
        [min(loc(:,1)) (min(loc(:,1))+ max(loc(:,1)))/2 max(loc(:,1)) xopt]
