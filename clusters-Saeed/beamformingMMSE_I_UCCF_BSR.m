%Sum-rate with BASR AP selection
function [Cmmse1, A_APs_selected]=beamformingMMSE_I_UCCF_BSR(channelCellfree,channelCellfree_tilde,Es,sigma,K)
% Pmmse MMSE precoder
% Pzf ZF precoder
% Cmmse MMSE sum-rate
% Czf ZF sum-rate
[a b]=size(channelCellfree*channelCellfree');
I=eye(a);
channelCellfree_hat=channelCellfree-channelCellfree_tilde;

%MMSE
f=sqrt(trace((((channelCellfree*channelCellfree')+sigma^2*I))^-2*(channelCellfree*channelCellfree'))/Es);
P1mmse=1/f*channelCellfree'*pinv(channelCellfree*channelCellfree'+sigma^2*I);

A_APs_selected=APselection_BSR(channelCellfree,channelCellfree_tilde,Es,sigma);
for ii=1:K
    P1mmse_a(:,ii)=A_APs_selected(:,:,ii)*P1mmse(:,ii);
end


Interf1=sigma^2*I+channelCellfree_tilde*P1mmse_a*P1mmse_a'*channelCellfree_tilde';
Cmmse1=log2(det(I+(channelCellfree_hat*P1mmse_a*P1mmse_a'*channelCellfree_hat')*pinv(Interf1)));









