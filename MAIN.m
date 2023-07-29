clear all;
close all;
K=8;M=64;
N=K*M;  
lambda=3e8/100e9;
d=lambda/2;
r=100
pos_DUA=[0:N-1]*d;

P=2;Q=5;
S1_CMS=[];S2_CMS=[];
for q=0:Q-1
    S1_CMS=[S1_CMS q*P*M*d+[0:M-1]*d] 
end
for p=0:2*P-1
    S2_CMS=[S2_CMS p*Q*M*d+[0:M-1]*d] 
end
S_CMS=[S1_CMS,S2_CMS(M+1:end)]; 
S_CMS=sort(S_CMS,'ascend');

S1_NMS=[];S2_NMS=[];
for q=0:K/2-1
     S1_NMS=[S1_NMS q*M*d+[0:M-1]*d];
     S2_NMS=[S2_NMS K/2*M*d+q*M*d*(K/2+1)+[0:M-1]*d];
end
S_NMS=[S1_NMS,S2_NMS]; 
S_NMS=sort(S_NMS,'ascend');

NR=[0,1,4,9,15,22,32,34];
S_NRMS=[];
for k=1:K
    S_NRMS=[S_NRMS NR(k)*M*d+[0:M-1]*d];
end

% DoF_DUA=DoF(pos_DUA);
% DoF_CMS=DoF(S_CMS);
% DoF_NMS=DoF(S_NMS);
% DoF_NRMS=DoF(S_NRMS);



dS_DUA=Channel_rank(N,pos_DUA,r,lambda);
dS_CMS=Channel_rank(N,S_CMS,r,lambda);
dS_NMS=Channel_rank(N,S_NMS,r,lambda);
dS_NRMS=Channel_rank(N,S_NRMS,r,lambda);

 co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
 
figure
semilogy(dS_DUA,'sk-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 7.2)
hold on
semilogy(dS_CMS,'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.8)
hold on
semilogy(dS_NMS,'dk-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.8)
hold on
semilogy(dS_NRMS,'ok-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5) 
grid on

 axis([1,30,1e-3,1])
 
 lgh=legend('DUA','Proposed CMS','Proposed NMS','Proposed NRMS')
set(lgh,'interpreter','latex','fontsize',14);
xlabel('Channel Rank','interpreter','latex','fontsize',14)
ylabel('Singular Value','interpreter','latex','fontsize',14)
b=[];b_theta=[];b_r=[];
theta=pi/3; 
for i=1:K*M
    nt=pos_DUA(i);
rnt=r^2-2*nt*r*sin(theta)+nt^2;
b(i)=1/sqrt(K*M)*exp(-1i*2*pi/lambda*sqrt(rnt));
b_theta(i)=1i*2*pi/lambda/sqrt(K*M)*exp(-1i*2*pi/lambda*sqrt(rnt))*nt*r*cos(theta)/sqrt(rnt);
b_r(i)=1i*2*pi/lambda/sqrt(K*M)*exp(-1i*2*pi/lambda*sqrt(rnt))*(nt*sin(theta)-r)/sqrt(rnt);
end

for j=1:40
    r=j*3
[CRB_theta_DUA,CRB_r_DUA(j)] = CRB(N,lambda,theta,r,pos_DUA);
[CRB_theta_CMS,CRB_r_CMS(j)] = CRB(N,lambda,theta,r,S_CMS);
[CRB_theta_NMS,CRB_r_NMS(j)] = CRB(N,lambda,theta,r,S_NMS);
[CRB_theta_NRMS,CRB_r_NRMS(j)] = CRB(N,lambda,theta,r,S_NRMS);
end 
figure
semilogy(3:3:120,CRB_r_DUA,'sk-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 7.2)
hold on
semilogy(3:3:120,CRB_r_CMS,'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.8)
hold on
semilogy(3:3:120,CRB_r_NMS,'dk-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.8)
hold on
semilogy(3:3:120,CRB_r_NRMS,'ok-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5) 
grid on

 axis([3,120,1e-5,1])
  lgh=legend('DUA','Proposed CMS','Proposed NMS','Proposed NRMS')
set(lgh,'interpreter','latex','fontsize',14);
xlabel('Range $$r$$ [meters]','interpreter','latex','fontsize',14)
ylabel('Root Range CRB','interpreter','latex','fontsize',14)