function [CRB_theta,CRB_r] = CRB(N,lambda,theta,r,array)

S_theta2=0;S_theta=0;S_r=0;S_r2=0;S_theta_r=0;
for i=1:N
       nt=array(i);
S_theta2=S_theta2+nt^2/(r^2-2*r*nt*sin(theta)+nt^2);
S_theta=S_theta+nt/sqrt(r^2-2*nt*r*sin(theta)+nt^2);
S_r=S_r+(nt*sin(theta)-r)/sqrt(r^2-2*nt*r*sin(theta)+nt^2);
S_r2=S_r2+(r^2-2*nt*r*sin(theta)+nt^2*sin(theta)^2)/(r^2-2*r*nt*sin(theta)+nt^2);
S_theta_r=S_theta_r+(nt*(nt*sin(theta)-r))/(r^2-2*r*nt*sin(theta)+nt^2);
end
 
chi=4*pi^2*r^2*cos(theta)^2/(lambda^2);

Q_matrix=[chi*(S_theta2/N-S_theta^2/((N)^2)), chi/r/cos(theta)*(S_theta_r/N-S_theta*S_r/((N)^2));...
         chi/r/cos(theta)*(S_theta_r/N-S_theta*S_r/(N^2)), chi*(S_r2/N-S_r^2/((N)^2))];
CRB_theta=sqrt(1/2/N*Q_matrix(2,2)/det(Q_matrix)) ;
CRB_r=sqrt(1/2/N*Q_matrix(1,1)/det(Q_matrix)) ;
end

