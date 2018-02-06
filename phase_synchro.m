%Phase Synchronization Main function
function [ ffmd_sync ] = phase_synchro( X,Ws,Lp, ch1id,ch2id )

    Lw=1;
    Z=Ws;
    X1=X;
    data = exp(1j * X1);

    windows=ceil((length(X1)-Ws+1)/Lp); 
    summed_cmpx_val=0;

    for i=1:windows
        fprintf('T_Windows: %d_%d of Channel Pairs %d_%d \n',i,windows,ch1id,ch2id);
        X=data(Lw:Z);
        Lw=Lw+Lp;
        Z=Z+Lp;
        summed_cmpx_val = sum(X);
        sync(i) = summed_cmpx_val;
    end

    ffmd_sync = (abs(sync)/Ws).^2; %[0 1]

end


