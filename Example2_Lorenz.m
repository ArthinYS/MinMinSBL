%% Example2: Learning the chaotic Lorenz system

clear; clc;close all;

addpath('Tools')

%% Generate Data
polyorder = 4;usesine = 1;n = 3;

% Lorenz's parameters (chaotic)
sigma = 10;  beta = 8/3; rho = 28;

% Initial condition
x0=[-8; 7; 27];  

% Integrate
dt = 0.001;tspan=[dt:dt:65];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);

%% Compute derivative

for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),sigma,beta,rho);
end

sp1=0;sp2=0;sp3=0;sp4=0;sp5=0; % record the number of successful trials

%% Add Gaussian noise to derivates

for j=1:100

 if mod(j,10)==0 
        fprintf('This is the %d-th iteration \n',j) 
 end



nls=0.001;
dx=dx+nls*randn(size(dx,1),size(dx,2));


subsam=1000;
sdx=dx(1:subsam:end,:);


%% Pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
sTheta=Theta(1:subsam:end,:);
m = size(Theta,2);
WT=zeros(95,3);
WT(2,1)=1;WT(3,1)=1;
WT(2,2)=1;WT(3,2)=1;WT(7,2)=1;
WT(4,3)=1;WT(6,3)=1;
Iter_Max=500;

%%  Our_SBL

for i=1:3
    [w,lam1,Iter_Num1]=Our_SBL(sTheta,sdx(:,i),Iter_Max);
    w1(:,i)=w;
end

ourw=w1;ourw(ourw~=0)=1;

if ourw-WT==0
    sp1=sp1+1;
end

%% Mackay_SBL

for i=1:3
   [w,lam2,Iter_Num2]=Mackay_SBL(sTheta,sdx(:,i),Iter_Max);
   w2(:,i)=w;
end


Mackayw=w2;Mackayw(Mackayw~=0)=1;

if Mackayw-WT==0
    sp2=sp2+1;
end


%% EM_SBL

for i=1:3
   [w,lam3,Iter_Num3]=EM_SBL(sTheta,sdx(:,i),Iter_Max);
    w3(:,i)=w;
end


EMw=w3;EMw(EMw~=0)=1;

if EMw-WT==0
    sp3=sp3+1;
end


%% IR_SBL

for i=1:3
    [w,Iter_Num4]=IR_SBL(sTheta,sdx(:,i),lam1,Iter_Max);
     w4(:,i)=w;
end

IRw=w4;IRw(IRw~=0)=1;

if IRw-WT==0
    sp4=sp4+1;
end



%% VI_SBL

for i=1:3
     [w,lam5]=VI_SBL(sTheta,sdx(:,i),1000);
      w5(:,i)=w;
end


VIw=w5;VIw(VIw~=0)=1;

if VIw-WT==0
    sp5=sp5+1;
end

end
%% Basis functions

 dictionary=poolDataLIST({'x','y','z'},w,n,polyorder,usesine)




