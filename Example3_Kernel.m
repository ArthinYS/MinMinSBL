%% Example 3: Sparse kernel regression 

clear;close all;clc
addpath('Tools')

%% Data

TN=1000;  Iter_Max=500;
Data=xlsread('winequality.xlsx');


%% Data preprocessing

Input=Data(:,1:end-1);

   for i=1:size(Input,2)
    Input(:,i)=2*(Input(:,i)-min(Input(:,i)))/(max(Input(:,i))-min(Input(:,i)))-1;
   end 

Output=Data(:,end);
Train_Input=Input(1:TN,:);Train_Output=Output(1:TN,:);
Test_Input=Input(TN+1:end,:);Test_Output=Output(TN+1:end,:);

%% Different kernels 
% 1: Exponential Kernel  2: Martern-3/2 kernel 
% 3: Linear kernel       4: Gaussian kernel


Ker_Ind=1;

%% Generate the dictionary using kernel functions

Phi(:,1)=ones(TN,1);
for i=1:size(Train_Input,1)
   for j=1:size(Train_Input,1) 
       
   
     if Ker_Ind==1
       % Exponential Kernel
         Phi(i,j+1)=exp(-norm(Train_Input(i,:)-Train_Input(j,:),2));
     
     elseif Ker_Ind==2
           % Martern-3/2 kernel
           nm=norm(Train_Input(i,:)-Train_Input(j,:),2);
           Phi(i,j+1)=(1+sqrt(3)*nm)*exp(-sqrt(3)*nm);

    elseif Ker_Ind==3
           % Linear kernel
           Phi(i,j+1)=Train_Input(i,:)*Train_Input(j,:)';
     else
          % Gaussian kernel
          nm=norm(Train_Input(i,:)-Train_Input(j,:),2)^2;
          Phi(i,j+1)=exp(-nm);

     end

    end
end

y=Train_Output;

%% Our_SBL
tic
[w1,lam1,Iter_Num1]=Our_SBL(Phi,y,Iter_Max);
toc
N1=size(nonzeros(w1),1);

%%  Mackay_SBL
 tic
[w2,lam2,Iter_Num2]=Mackay_SBL(Phi,y,Iter_Max);
 toc
N2=size(nonzeros(w2),1);

    %% EM_SBL
 tic
 [w3,lam3,Iter_Num3]=EM_SBL(Phi,y,Iter_Max);
 toc

N3=size(nonzeros(w3),1);

 %% IR_SBL
 tic
 [w4,Iter_Num4]=IR_SBL(Phi,y,lam1,Iter_Max);
 toc
 N4=size(nonzeros(w4),1);

 %% VI_SBL
 tic
 [w5,lam5]=VI_SBL(Phi,y,1000);
 toc
 N5=size(nonzeros(w5),1);

 %% Compute the MRE of each SBL algorithm on testing dataset

TNT=size(Test_Output,1);
PhiP(:,1)=ones(TNT,1);
for i=1:size(Test_Input,1)
    for j=1:size(Train_Input,1)
  
  
     if Ker_Ind==1
          % Exponential Kernel
          PhiP(i,j+1)=exp(-norm(Test_Input(i,:)-Train_Input(j,:),2));

      elseif Ker_Ind==2
           % Martern-3/2 kernel  
           nm=norm(Test_Input(i,:)-Train_Input(j,:),2);
           PhiP(i,j+1)=(1+sqrt(3)*nm)*exp(-sqrt(3)*nm);

    elseif Ker_Ind==3
           % Linear kernel
           PhiP(i,j+1)=Test_Input(i,:)*Train_Input(j,:)';
    else
           % Gaussian kernrl
           nm=norm(Test_Input(i,:)-Train_Input(j,:),2)^2;
           PhiP(i,j+1)=exp(-nm);
     end
    end
end

MRE1=sum(abs((Test_Output-PhiP*w1)./Test_Output))/TNT;
MRE2=sum(abs((Test_Output-PhiP*w2)./Test_Output))/TNT;
MRE3=sum(abs((Test_Output-PhiP*w3)./Test_Output))/TNT;
MRE4=sum(abs((Test_Output-PhiP*w4)./Test_Output))/TNT;
MRE5=sum(abs((Test_Output-PhiP*w5)./Test_Output))/TNT;

TN1=TN+1;

%% compute d=sqrt(MRE^2+Sparsity^2)

d1=sqrt(MRE1^2+(N1/TN1)^2);
d2=sqrt(MRE2^2+(N2/TN1)^2);
d3=sqrt(MRE3^2+(N3/TN1)^2);
d4=sqrt(MRE4^2+(N4/TN1)^2);
d5=sqrt(MRE5^2+(N5/TN1)^2);


