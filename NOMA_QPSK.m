clc;
clear all;
close all;
N=1000;
snr=15;

 x1=randi([0,1],N,1); %bits generated user1
 x2=randi([0,1],N,1);%bits generted user2

QPSKmod = comm.QPSKModulator('BitInput',true); 
QPSKdemod = comm.QPSKDemodulator('BitOutput',true); 
qpsk_const=reshape(constellation(QPSKmod),1,4)

x1_modulated=reshape(QPSKmod(x1),[1,N/2]) ;%modulated symbols user1
x2_modulated=reshape(QPSKmod(x2),[1,N/2]) ;%modulated symbols user2

h1=(1/sqrt(2))*(randn(1,N/2)+1i*randn(1,N/2));%channel coefficient u1
h2=(1/sqrt(2))*(randn(1,N/2)+1i*randn(1,N/2));%channel coefficient u2

y=h1.*x1_modulated+h2.*x2_modulated ;%received symbol

[z,var] = awgn(y,snr,'measured'); %awgn noise addition


%successive interference cancellation
decoded_symbols_u1=zeros(1,N/2);
decoded_symbols_u2=zeros(1,N/2);
decoded_bits_u1=zeros(1,N);
decoded_bits_u2=zeros(1,N);
for k=1:N/2
if(abs(h1(k))>abs(h2(k)))
    %decode user 1 symbol first
decoded_symbols_u1(k)=ml_decoder(qpsk_const,h1(k),z(k)); %directly decode user1 symbols
z_u2(k)=z(k)-h1(k)*decoded_symbols_u1(k); %interference cancellation
decoded_symbols_u2(k)=ml_decoder(qpsk_const,h2(k),z_u2(k));%decode user2 symbols
else
 decoded_symbols_u2(k)=ml_decoder(qpsk_const,h2(k),z(k));%decode user2 symbols
z_u1(k)=z(k)-h2(k)*decoded_symbols_u2(k);%interference cancellation
decoded_symbols_u1(k)=ml_decoder(qpsk_const,h1(k),z_u1(k)); %directly decode user1 symbols

end


end
decoded_symbols_u1;
decoded_symbols_u2;
decoded_bits_u1=QPSKdemod(reshape(decoded_symbols_u1,[N/2,1]));  %decoded bits user1
decoded_bits_u2=QPSKdemod(reshape(decoded_symbols_u2,[N/2,1])); %decoded bits user2
ber_u1=biterr(x1,decoded_bits_u1)/N %ber user1
ber_u2=biterr(x2,decoded_bits_u2)/N %ber user2


%joint decoding
for l=1:N/2
    [u1(l),u2(l)]=joint_ml_decode(qpsk_const,z(l),h1(l),h2(l));
    
end
u1;
u2;
bits_u1=QPSKdemod(reshape(u1,[N/2,1]));
bits_u2=QPSKdemod(reshape(u2,[N/2,1]));
    
bits_u1;
bits_u2;
ber_joint_u1=biterr(x1,bits_u1)/N
ber_joint_u2=biterr(x2,bits_u2)/N 


function y=ml_decoder(const,h,rx_symbols)
    const_user=h*const;
distance=zeros(1,length(const));
for i=1:length(const)
    distance(i)=(rx_symbols-const_user(i))*conj(rx_symbols-const_user(i));
end
min_distance=min(distance);
index=(distance==min_distance);
y=const(index);
end

function [y1,y2]=joint_ml_decode(const,z,h1,h2)
const_u1=h1*const;
const_u2=h2*const;
k=length(const);
dist=zeros(k);
for i=1:length(const)
    for j=1:length(const)
       dist(i,j)=(const_u1(i)+const_u2(j)-z)*conj(const_u1(i)+const_u2(j)-z);
    end
end
min_dist=min(min(dist));
[r,c]=find(dist==min_dist);
y1=const(r);
y2=const(c);
end

