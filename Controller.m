function [e,r,Tau, theta_est_dot, theta_tilda] = controller(t,q,q_dot,theta_est)
%Raghavv Goel | 2016179
c2=cos(q(2));
s2=sin(q(2));
alpha=10;
k=10;


p1=3.473;
p2=0.196;
p3=0.242;
fd1=5.3;
fd2=1.1;

theta = [p1;p2;p3;fd1;fd2];

qd=[0.5*sin(t);2*cos(t/4)];
qd_dot=[0.5*cos(t);-0.5*sin(t/4)];
qd_ddot=[-0.5*sin(t);-0.25*cos(t/4)];
e=q-qd;     %Tracking Error
e_dot=q_dot-qd_dot;

r=e_dot+alpha*e;         %Filtered Tracking Error

y = 10*eye(5);
%----------------------------------------------------------%
%M=[p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];% Inertia matrix
%Vm=[-p3*s2*q_dot(2) -p3*s2*(q_dot(1)+q_dot(2));p3*s2*q_dot(1) 0];%Centripetal coriolis matrix
%fd=[fd1 0;0 fd2]; % Friction matrix

z1 = [-qd_ddot(1) + alpha*e_dot(1) ; 0];
z2 = [-qd_ddot(2) + alpha*e_dot(2) ; -qd_ddot(1) -  qd_ddot(2) + alpha*(e_dot(1) + e_dot(2))];
z3 = [2*c2*(-qd_ddot(1) + alpha*e_dot(1)) + c2*(-qd_ddot(2) + alpha*e_dot(2)) + s2*q_dot(2)*(q_dot(1) + r(1)) + s2*(q_dot(1) + q_dot(2))*(q_dot(2) + r(2)) ;
      c2*(-qd_ddot(1) + alpha*e_dot(1)) - s2*q(1)*(q_dot(1) + r(1))  ];
z4 = [-q_dot(1); 0];
z5 = [0; -q_dot(2)];

Z = [z1 z2 z3 z4 z5];
%Z = [-qd_ddot(1) + alpha*q_dot(1) - alpha*qd_dot(1), -qd_ddot(2) + alpha*(q_dot(2) - qd_dot(2)), -2*c2*qd_ddot(1) - c2*qd_ddot(2) + 2*s2*q_dot(1)*q_dot(2) + s2*q_dot(2)*q_dot(2) + 2*alpha*c2*(q_dot(1) - qd_dot(1)) + alpha*c2*(q_dot(2) - qd_dot(2)), q_dot(1) , 0 ; 0, -qd_ddot(1) - qd_ddot(2) + alpha*(q_dot(1) + q_dot(2) - qd_dot(1) - qd_dot(2)),  -c2*qd_ddot(1) - s2*q_dot(1)*q_dot(1) + alpha*(q_dot(1) - qd_dot(1)), 0, q_dot(2)];

theta_est_dot = y*Z'*r;

%theta_est = integral(theta_est_dot,0,inf);

Tau = -Z*theta_est - k*r - e;
%tau=-k*r+Vm*q_dot+fd*q_dot+M*qd_ddot-alpha*M*e_dot-Vm*r-e;         % Controller
theta_tilda = theta - theta_est;

%figure
%plot(theta_tilda, t);
