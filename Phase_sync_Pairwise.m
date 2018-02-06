% Frobenius norm PHASE COHERENCE
% CODED BY: PUNEET DHEER(pd)
% Date: 05-09-2017
% INPUT:
% channel_matrix= input data row wise channel
% Ws = Window size (in sample point)
% Lp = shift window by sample point 
% 
% OUTPUT:
% ps_matrix = pair wise phase coherence based on Frobenius norm 
% ensemble_measure = common synchronization among all the channels 
%%
function [ ps_matrix,ensemble_measure] = Phase_sync_Pairwise( channel_matrix,Ws,Lp )
% outputfile = 'result.png';
tic
n=size(channel_matrix,1);%no. of channels
% windowed_eigen_vector=[];
ps_matrix=[];
ps_matrix_across=zeros(n,n);

numChannels = size(channel_matrix, 1);

for i = 1:numChannels
    Phase1 = angle(hilbert(channel_matrix(i, :)));
    Phase(i,:) = mod(Phase1,2*pi);
    %Phase(i, :) = atan2(imag(hilbert(channel_matrix(i, :))),real(hilbert(channel_matrix(i, :))));
end


for i=1:n
    for j=i:n
        i
        j
        %tic
        phase_diff=mod((Phase(i,:)-Phase(j,:)),2*pi);
        [P_S]= phase_synchro(phase_diff,Ws,Lp,i,j);
        
        %toc
        ps_matrix=[ps_matrix;P_S]; %storing pairwise 
        
        
    end
end

p=size(ps_matrix,2); %no. of windows


for i=1:p %no. of windows
    a=0;
    for j=1:n
        for k=j:n
            a=a+1;
            ps_matrix_across(j,k)=ps_matrix(a,i);     
            ps_matrix_across(k,j)=ps_matrix(a,i);
        end
    end
    diagonal=diag(ps_matrix_across); %just to check
%    %Note: sum of eigen value should be equal to the no. of Channels.
%    %sum of diagonal of windowed PLV Matrix is also equal to no. of Channels.
%    windowed_matrix=ps_matrix_across;
%    V=sort(eig(ps_matrix_across));
%     sum(V)
%    windowed_eigen_vector=[windowed_eigen_vector V];
%    norm_eig = V./sum(V);
%    S_esti(i)=real(1+(sum(norm_eig.*log(norm_eig))/log(n)));
%     
   fro_norm(i) = norm(ps_matrix_across ,'fro');
     
   
   ps_matrix_across=zeros(n,n);
  
    
end
ensemble_measure = sqrt((fro_norm.^2 - n)./(n^2 - n));

toc




end

