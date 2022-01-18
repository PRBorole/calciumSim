function fval = temp(t,x)

A = x(1);
B = x(2);
C = x(3);
D = x(4);


kf = 0.2;
kr = 0.1;
% dt = 0.00001;
% 
% A_(1) = A;
% B_(1) = B;
% C_(1) = C;
% D_(1) = D;

% for i = 2:10^3
    dA = -kf*A*B + kr*C*D;
%     A_(i) = A + dt*dA;
    
    dB = -kf*A*B + kr*C*D;
%     B_(i) = B + dt*dB;
    
    dC = kf*A*B - kr*C*D;
%     C_(i) = C + dt*dC;
    
    dD = kf*A*B - kr*C*D;
%     D_(i) = D + dt*dD;
%     
%     A = A + dt*dA;
%     B = B + dt*dB;
%     C = C + dt*dC;
%     D = D + dt*dD;
% end
% 
% figure;
% plot(1:1000, A_);
% figure;
% plot(1:1000, C_);



fval = [dA;dB;dC;dD];




