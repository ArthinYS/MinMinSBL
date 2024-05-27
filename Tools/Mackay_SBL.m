function [w,lam,Iter_Num]=Mackay_SBL(Phi,Y,Iter_Max)

%% Tipping_SBL1 with MacKay updates

%% Initializations 
[N,M]=size(Phi); ind=1:M;
threshold1=1e-6; threshold2=1e-4; 
Gamm=10*ones(M,1); lam=0.1;

Iter_Num=0;
%% Learning loop
N1=M;
for i=1:1:Iter_Max
    Iter_Num=Iter_Num+1;
    Gamm_old=Gamm;
    Sigm_1in =  Phi'*Phi+diag(lam*Gamm);
   % Sigm_2=lam*eye(N)+Phi*diag(Gamm)*Phi';
    mu=Sigm_1in\(Phi')*Y;
    Sigma_1=lam*(Sigm_1in\eye(N1));
    ga=diag(eye(N1))-Gamm.*diag(Sigma_1);
    Gamm=ga./(mu.^2);
    lam=norm(Y-Phi*mu,2)^2;lam=lam/(N-sum(ga));
    if max(max(abs(1./Gamm-1./Gamm_old)))<threshold1
        break;
    end
    ind_dele=find(1./Gamm<threshold2);
    Gamm(ind_dele)=[];
    Phi(:,ind_dele)=[]; 
    N1=size(Phi,2);
    ind(ind_dele)=[];
end

 Sigm_1in =  Phi'*Phi+diag(lam*Gamm);x_opt=Sigm_1in\(Phi')*Y;
 w=zeros(M,1);w(ind)=x_opt;
end