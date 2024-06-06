
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% This example demonstrates that all the employed SBL codes are indeed 
%% effective in finding sparse solutions with appropriate experimental settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  It should be highly noted that 
%% you have to download the CVX toolbox from https://cvxr.com/ to run IR_SBL%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Here, we use all the SBL algorithms to identify y=x^3 from the data 
%% given the basis functions: 1,x,x^2,x^3,x^4,x^5.

clear;clc;close all

%% Generate data
x=-20:1:20; x=x'; 
y=x.^3+0.001*randn;

%% Basis function
fx=[ones(size(x,1),1),x,x.^2,x.^3,x.^4,x.^5];
fx_name={'1','x','x.^2','x^3','x^4','x^5'};
Iter_Max=500;

%% Our_SBL
  [w1,lam1,Iter_Num1]=Our_SBL(fx,y,Iter_Max);
%% Mackay_SBL
  [w2,lam2,Iter_Num2]=Mackay_SBL(fx,y,Iter_Max);
%% EM_SBL
  [w3,lam3,Iter_Num3]=EM_SBL(fx,y,Iter_Max);
%% IR_SBL
  [w4,Iter_Num4]=IR_SBL(fx,y,lam1,Iter_Max);
%%  VI_SBL
  [w5]=VI_SBL(fx,y,1000);
