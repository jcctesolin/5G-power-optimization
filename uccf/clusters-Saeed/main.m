



clc
clear all
close all


% ----- Modulation scheme is 4-QAM
M=2;k_QPSK=sqrt(2);a=sqrt(2);delta=3/4;

% Linear MMSE and ZF equalizer

% Snr in dB
EbN0=[0:5:25];



%Set the side length of the simulation area
squareLength = 400; %for cell free



% Number of APs in cell-free network
nbrOfAPs = 16;

%Number of antennas per AP
N = 4;

% Trial number
T=100;

%Number of realizations of the random UE locations
nbrOfSetups =100;


%Number of UEs in the simulation setup
K = 32;



% Symbol energy 
Ps = 1;

% Transmitted power
Es=K*Ps;
I=eye(K);





% cell-free without scheduling Imperfect CSI
   C_MMSE_IP1=zeros(size(EbN0));%without power allocation

    
%--------------------------------------------------------------------


%-------------------------UCCF without user scheduling for different AP selection criteria

C_MMSE_a1=zeros(size(EbN0));%without power allocation 



C_MMSE_basr1=zeros(size(EbN0));%without power allocation 


%Set the AP locations for the small-cell and cell-free setups
APperdim = sqrt(nbrOfAPs);
APcellfree = linspace(squareLength/APperdim,squareLength,APperdim)-squareLength/APperdim/2;
APcellfree = repmat(APcellfree,[APperdim 1]) + 1i*repmat(APcellfree,[APperdim 1])';



ro_f=1;

% Counting number
ns0=0;
% Starting Programm
for EN=EbN0
    %The programming indicator
    ns0=ns0+1
    for trial=1:T
        for pack = 1:nbrOfSetups
            
UElocations =(rand(nbrOfSetups,K)+1i*rand(nbrOfSetups,K))*squareLength;


            %Generate the channel matrix in the cell-free setup
            [channelCellfree_I, ~,channelCellfree_tilde]=channelcellfree_N(APcellfree,UElocations,nbrOfAPs,N,K,pack);
            

            
          
EbN0_cnv=exp(EN*log(10)/10);



          
           
            
           


                                sig = sqrt((trace(channelCellfree_I*channelCellfree_I'))/(EbN0_cnv*K));

                            
        
        
%------------------------------------------------------------------------------------------------------
            % sum rates with ZF and MMSE precodings for CoMP and CF
            % networks without scheduling for perfect CSI and Imperfect CSI

           

            
  [C_mmse_IP1, C_mmse_IP1_EPL]=beamformingMMSE_I(channelCellfree_I,channelCellfree_tilde,Es,sig);
            C_MMSE_IP1(ns0)=C_MMSE_IP1(ns0)+C_mmse_IP1;



%--------------------------------------------------------------------------------------------------------

 %------------------UCCF simulation MMSE precoder only----------------
[~, largeScaleFading, ~]=channelcellfree_N(APcellfree,UElocations,nbrOfAPs,N,K,pack);
            
           
                              [Cmmse1_a]=beamformingMMSE_I_UCCF(channelCellfree_I,channelCellfree_tilde,Es,sig,K,largeScaleFading,N);

            C_MMSE_a1(ns0)=C_MMSE_a1(ns0)+Cmmse1_a;%without power allocation

           
            %------------------end of UCCF simulation MMSE precoder only----------------


            

             %------------------UCCF simulation MMSE precoder only, Boosted acheivable rate criterion----------------


            [Cmmse1_basr]=beamformingMMSE_I_UCCF_BSR(channelCellfree_I,channelCellfree_tilde,Es,sig,K);

            C_MMSE_basr1(ns0)=C_MMSE_basr1(ns0)+Cmmse1_basr;%without power allocation
            %------------------end of UCCF simulation MMSE precoder only, Boosted acheivable rate criterion----------------

              

            
          
            
        end
    end
    
end



% cell-free without scheduling Imperfect CSI
C_MMSE_IP1=C_MMSE_IP1/(nbrOfSetups*T);
    



%-------------------------UCCF without user scheduling for different AP selection criteria

C_MMSE_a1=C_MMSE_a1/(nbrOfSetups*T);%without power allocation 

C_MMSE_basr1=C_MMSE_basr1/(nbrOfSetups*T);%without power allocation 





%--------------------
%--------------------------
figure;
semilogy(EbN0, C_MMSE_a1,'s-b','LineWidth',2,'Markersize',8);
hold on;
semilogy(EbN0,C_MMSE_basr1,'o-k','LineWidth',2,'Markersize',8);
semilogy(EbN0,C_MMSE_IP1,'^-r','LineWidth',2,'Markersize',8);




xlabel('SNR-dB');
ylabel('Sum-Rate');
title('Sum Rate wituout rsource allocation, Comparison of LSF fading and BSR criteria');
legend('UCCF-LSF','UCCF-BSR','CF')
grid





