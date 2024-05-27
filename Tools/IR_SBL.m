function [w,Iter_Num]= IR_SBL(Phi,Y,lambda,Iter_Max)

%% Initializations 

[M,N]=size(Phi);  ind=1:N;

threshold1=1e-6; threshold2=1e-4; 

Gamm=0.1*ones(N,1);
Iter_Num=0;
% initialisation of the variables
U=ones(N,1);

%fprintf(1, 'Finding a sparse feasible point using l1-norm heuristic ...\n');

for iter=1:Iter_Max
    N1=size(Phi,2);
    Iter_Num=Iter_Num+1;
    Gamm_old=Gamm;
 %   fprintf('This is round %d \n', iter);
% cvx_clear
    cvx_begin quiet
    cvx_solver sedumi   %sdpt3
    variable W(N1)
    minimize    (lambda*norm( U.*W, 1 )+ 0.5*sum((Phi* W-Y).^2) )
    %                 subject to
    %                           W.^2-ones(101,1)<=0;
    cvx_end
    
    Gamm=(U.^-1).*W;

    if max(max(abs(Gamm-Gamm_old)))<threshold1
        break;
    end

    ind_dele=find(Gamm<threshold2);

    Gamm(ind_dele)=[];
    W(ind_dele)=[];
    Phi(:,ind_dele)=[]; 

    Usss=size(Phi,1);
    if Usss==0
        break;
    end

    N1=size(Phi,2);
    ind(ind_dele)=[];

    Dic0=lambda*eye(M)+Phi*diag(Gamm)*Phi';
    %%
    UU=diag(Phi'*(Dic0\Phi));
    U=abs(sqrt(UU));
    Us=size(U,1);
    if Us==0
        break;
    end
end
 w=zeros(N,1);
 w(ind)=W;   
end

