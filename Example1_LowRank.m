
%% Example 1: Sparse signal recovery

%%%% We conduct 100 simulation trials with the SNR values 
%%%% ranging from 0 to 35 dB in steps of 5 dB.


%%  
clear;close all;clc

addpath('Tools')

sp1=0;sp2=0;sp3=0;sp4=0;sp5=0; % record the number of successful trials

tic

for j=1:100

    if mod(j,10)==0 
        fprintf('This is the %d-th iteration \n',j) 
    end

%% SNR

 snr=5;

%% Experimental setting

n=60; m=100; ms=4; Iter_Max=100; 

w = zeros(m,1);
ids = randperm(m,ms);

%% Record the nonzero terms

Ind=zeros(1,m)';

for i=1:ms
    Ind(ids(i))=1;
end

%% Generate sparse signal
% ss=1: spike signal; ss=2: Gaussian signal; ss=3: Uniformly distributed signal

ss=1;
if ss==1
   val = randi(2,ms,1)-1.5 ;
   w(ids) = val*2;

elseif ss==2
    val=randn(ms,1);
    w(ids) = val;

else
    val=2*rand(ms,1)-1;
    w(ids) = val;

end

%% Generate low-rank matrices

r = 40; % the rank of the matrix 

Phi_in = rand(n,m); % initial matrix

[U, S, V] = svd(Phi_in); % SVD decomposition

S_truncated = S; 
S_truncated(r+1:end, r+1:end) = 0; % truncate singular values

Phi = U * S_truncated * V'; % low-rank matrix


%% Nosiy output

   y = Phi * w  ;

   y = awgn(y, snr, 'measured');

%% Our_SBL

[w1,lam1,Iter_Num1]=Our_SBL(Phi,y,Iter_Max);

% Successful probability
nw1=w1;nw1(nw1~=0)=1;
if sum(Ind-nw1)==0
    sp1=sp1+1;
end

 %% Mackay_SBL

 [w2,lam2,Iter_Num2]=Mackay_SBL(Phi,y,Iter_Max);

% Success probability
nw2=w2;nw2(nw2~=0)=1;
if sum(Ind-nw2)==0
    sp2=sp2+1;
end

 %% EM_SBL

 [w3,lam3,Iter_Num3]=EM_SBL(Phi,y,Iter_Max);

% Success probability
nw3=w3;nw3(nw3~=0)=1;
if sum(Ind-nw3)==0
    sp3=sp3+1;
end

 %% IR_SBL
 [w4,Iter_Num4]=IR_SBL(Phi,y,lam1,Iter_Max);

 % Success probability
   nw4=w4;nw4(nw4~=0)=1;
if sum(Ind-nw4)==0
    sp4=sp4+1;
end

  %% VI_SBL
 [w5,lam5]=VI_SBL(Phi,y,1000);
 
 % Success probability
 nw5=w5;nw5(nw5~=0)=1;
 if sum(Ind-nw5)==0
     sp5=sp5+1;
 end
end
toc

