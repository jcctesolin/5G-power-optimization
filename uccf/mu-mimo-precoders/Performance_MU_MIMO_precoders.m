clear all;
% Number of transmissions
T = 10;
snr = 0:5:25;
b_zf = zeros(T,length(snr));
b_mmse = zeros(T,length(snr));
for j = 1:T
   for i = 1:length(snr)
      [b_zf(j,i),b_mmse(j,i)] = MU_MIMO_DL_precoders(snr(i));
   end
end
bzf = mean(b_zf,1);
bmmse = mean(b_mmse,1);

semilogy(snr,bzf,'r-o',snr,bmmse,'b-s');
xlabel('SNR');
ylabel('BER');
legend('ZF-Linear','MMSE-Linear')
grid;