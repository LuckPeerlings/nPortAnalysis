function f = ScrambleFrequency(MEASUREMENTFREQUENCY,REF_CHAN)
ii = 1;
f_list = MEASUREMENTFREQUENCY;
while isempty(f_list) == 0
            N = length(f_list);        
            ii_rand = ceil(N*rand);%Generate random integer number
            f_rand(ii) = f_list(ii_rand);
            f_list(ii_rand) = []; 
            ii = ii+1;
end
f(:,1) = f_rand;


%The other frequencies are compared against the 
if size(REF_CHAN,2) == 2
    cc = 0; %Counter for the amount of tries to obtain an extra random vector
    ii = 1;
    f_list = MEASUREMENTFREQUENCY;
    while isempty(f_list) == 0
            N = length(f_list);
            ii_rand=ceil(N*rand);
            f_rand(ii) = f_list(ii_rand);
            if ( rem(f_rand(ii),f(ii)) == 0 || rem(f(ii),f_rand(ii)) == 0) ...
                || abs(f_rand(ii)/f(ii)-1) < 0.1
                disp('Frequency a multiple of each other or too close to each other')
                if N < 10
                    disp('Not enough frequencies left, starting a new iteration')
                    ii = 1;
                    f_list = MEASUREMENTFREQUENCY;
                end
                cc = cc + 1;
            else
                f_list(ii_rand) = [];
                ii = ii+1;
       end
       if cc > 100; error('Too many tries to find an appropriate frequency vector, check the settings'); end
       end
       f(:,2) = f_rand;
end
if size(REF_CHAN,2) > 2
    error('The script has not been coded to take care of more than two references')
end
