%Sum-rate without selection
function [Cmmse1, P1mmse]=beamformingMMSE_I(channelCellfree,channelCellfree_tilde,Es,sigma)
% Pmmse MMSE precoder
% Pzf ZF precoder
% Cmmse MMSE sum-rate
% Czf ZF sum-rate
channelCellfree_hat=channelCellfree-channelCellfree_tilde;
[a b]=size(channelCellfree*channelCellfree');
I=eye(a);
% sigma=sigma_CF(channelCellfree,EbN0_cnv,K);


%MMSE
f=sqrt(trace((((channelCellfree*channelCellfree')+sigma^2*I))^-2*(channelCellfree*channelCellfree'))/Es);
P1mmse=1/f*channelCellfree'*pinv(channelCellfree*channelCellfree'+sigma^2*I);

Interf1=sigma^2*I+channelCellfree_tilde*P1mmse*P1mmse'*channelCellfree_tilde';
Cmmse1=log2(det(I+(channelCellfree_hat*P1mmse*P1mmse'*channelCellfree_hat')*pinv(Interf1)));







