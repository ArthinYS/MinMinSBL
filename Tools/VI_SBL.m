function [w,sigma]=VI_SBL(Phi,Y,r)  %% r  can be set to be 1000       
         [N,M]=size(Phi);
         s=0;
         Ms=1:M;
         Mstar=M;
          while Mstar>1      
              Mstar=size(Phi,2);
              s=s+1;
              MS{s}=Ms;
             % fprintf('this is the %dth iteration\n',s)
              [thetak,EAA, LQ_final,sigma2,VkI]=VI(Y,Phi);
              ARD=1./diag(EAA);       
              logTARD=min(log(ARD))+(max(log(ARD))-min(log(ARD)))/r;
              indm=find(log(ARD)<logTARD);
              Phi(:,indm)=[];
              Ms(indm)=[];
              LQmax(s)=LQ_final;   
              Theta{s}=thetak;
              sigma_square{s}=sigma2;
          end
               [mannum, maxind]=max(LQmax);
               w=zeros(M,1);
               w(MS{1,maxind})=Theta{1,maxind};
              sigma=sqrt(sigma_square{maxind});
end 








%%
function [thetak,EAA,LQ_final,sigma2,VkI]=VI(Y,Phi)
%% Initialization
[N,M]=size(Phi);

if N<M
    nor=1e-1;
else
    nor=1e-4;
end

a0=1e-2;b0=1e-4;c0=1e-2;d0=1e-4;
ak=a0+N/2;bk=b0;ck=c0+1/2;dk=d0*ones(M,1);
%[u1,s1,v1]=svd(Phi'*Phi+nor*eye(M))
thetak=(Phi'*Phi+nor*eye(M))\Phi'*Y;
Vk=Phi'*Phi;
thr=1e-4;LQ(1)=-inf;
%%  main loop
     i=0;
  while 1
      i=i+1;
      EAA=diag(ck./(dk));
      ETT=(thetak.^2).*(ak/bk)+diag(Vk);
      VkI=Phi'*Phi+EAA;
      Vk=(VkI+nor*eye(M))\eye(M);
      Py=0;
      for j=1:N
         Py=Py+ Phi(j,:)'*Y(j);
      end
      thetak=Vk*Py;
      bk=b0+0.5*(sum(Y.^2)-thetak'*VkI*thetak);
      dk=d0+0.5*ETT;
      LQ(i+1)=0.5*log(det(Vk))-b0*(ak/bk)-0.5*(ak/bk*sum((Y-Phi*thetak).^2)+sum(diag(Phi*Vk*Phi')))-...
              ak*log(bk)-sum(ck*log(dk)); %% The constant term is not includeded herein.
          LQ_final=LQ(i+1);
      if   LQ(i+1)-LQ(i)<thr
          break;
      end
  end
  sigma2=bk/(ak-1);
end