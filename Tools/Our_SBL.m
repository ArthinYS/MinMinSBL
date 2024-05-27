function [w,lam,Iter_Num]=Our_SBL(Phi,Y,Iter_Max)

%% ============= Reference =============
% Yasen Wang, Junlin Li, Zuogong Yue, and Ye Yuan, 
% An iterative Min-Min optimization method for sparse Bayesian learning

%% ============= Author =============
%  Yasen Wang (arthinw@hust.edu.cn,arthinw@163.com)

%%%%%% Phi: Library matrix : N * M %%%%%%%%%%
% ===== N: The number of data points ====
% ===== M: The number of basis functions ===

%%%%   Y: Output vector : N*1    %%%%%%%

%% Initializations 
[N,M]=size(Phi); ind=1:M;
a0=1e-5; b0=1e-5;
threshold1=1e-6; threshold2=1e-4; 
Gamm=0.1*ones(M,1); lam=0.1;
 Iter_Num=0;
%% Learning loop
for i=1:1:Iter_Max
    Iter_Num=Iter_Num+1;
   % fprintf('This is the %d-th iteration\n',i)
    Gamm_old=Gamm;
    Sigm_1in =  Phi'*Phi+diag(lam*(1./Gamm));
    Sigm_2=lam*eye(N)+Phi*diag(Gamm)*Phi';
    z_opt1 =diag(Phi'*(Sigm_2\Phi))+(2*a0+2)*(1./Gamm);
    z_opt2=trace(Sigm_2\eye(N))+(2*a0+2)/lam;
    x_opt=(Sigm_1in\(Phi'))*Y;
    Gamm=sqrt((x_opt.^2+2*b0)./(z_opt1));
    lam=sqrt((norm(Y-Phi*x_opt,2)^2+2*b0)/z_opt2); 
    if max(max(abs(Gamm-Gamm_old)))<threshold1
        break;
    end
    ind_dele=find(Gamm<threshold2);
    Gamm(ind_dele)=[];
    Phi(:,ind_dele)=[];  
    ind(ind_dele)=[];
end


Sigm_1in =  Phi'*Phi+diag(lam*(1./Gamm));
x_opt=Sigm_1in\(Phi')*Y;
w=zeros(M,1);
w(ind)=x_opt;

end