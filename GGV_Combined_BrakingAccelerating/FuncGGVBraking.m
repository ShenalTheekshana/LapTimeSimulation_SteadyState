 function F = FuncGGVBraking(r,v,phi,sat)
%psi=0.1;
%phi=0.2;
g=9.81;  %gravitataional constant
ro=1.225;%air density
M=300;   %mass     
Rw=0.2325;   %wheel radius   
Iw=0.05;
gamma=5; %camber angle   
Af=1;    %frontal area
As=2;    %span area   
Cd=1;    %drag coefficient   
Cl=1.5;    %lift coefficient   
TmaxB=2000; %maximum brake torque  
TmaxM=400; %motor maximum torque  
vmax=30;

%Cd_front = 0.003228*(v^3) - 0.001536*(v^2) - 0.008723*v + 0.8459
%Cd_Rear = 0.0001991*(v^2) - 0.003648*v + 0.7156
%Cd_Under = 0.0001333*(v^3) - 0.009029*(v^2) + 0.1875*v - 0.594

%Cl_front = -0.02607*(v^2) + 0.07178*v + 3.641
%Cl_Rear = 0.00024*(v^2) + 0.00244*v + 3.016
%Cl_Under = -1.098*(v^2) + 0.887*v + 3.411

%A_front = 1.7731;
%A_rear = 1.6717;
%A_under = 1.7336;

%Cd = ((Cd_front*A_front)+(Cd_Rear*A_rear)+(Cd_Under*A_under))/(A_front + A_rear + A_under);
%Cl = ((Cl_front*A_front)+(Cl_Rear*A_rear)+(Cl_Under*A_under))/(A_front + A_rear + A_under);

%Cd = (((0.003228*(v^3) - 0.001536*(v^2) - 0.008723*v + 0.8459)*A_front)+((0.0001991*(v^2) - 0.003648*v + 0.7156)*A_rear)+((0.0001333*(v^3) - 0.009029*(v^2) + 0.1875*v - 0.594)*A_under))/(A_front + A_rear + A_under);
%Cl = ((( -0.02607*(v^2) + 0.07178*v + 3.641)*A_front)+((0.00024*(v^2) + 0.00244*v + 3.016)*A_rear)+((-1.098*(v^2) + 0.887*v + 3.411)*A_under))/(A_front + A_rear + A_under);


b = [5 1000];      %lateral tire force coefficients 
%    b1   b2
a = [3 1092 0.005];%longitudinal tire force coefficients
%    a1   a2   a15

%----------------------------------------------------------------------------------------
        
% r(1) r(2)  r(3) 
%  r   Fx    Fr  


 Fl=  0.5*Cl*ro*As*(v^2);                        %calculating the downforce (subscription of Fl means "lift" )
 Fz= (M*g+Fl)/4000;                              %calculating the total normal force (verticle force towards ground) in kN
 Fxmax=Fz*(b(1)*Fz + b(2));                      %calculating the maximum longitudinal force based on magic formula in N(D constant in magic formula)
 Fymax=Fz*(a(1)*Fz + a(2))*(1 - a(3)*(gamma^2)); %calculating the maximum lateral force based on magic formula in N(D constant in magic formula)
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %ay=r*cos(phi);
 %Fy=M*ay/4;  %Fy=M*(r*cos(phi))/4;  %calculating lateral tire force on one tire using F=ma
 %Fr= 0.01*Fz*1000;                  %calculating the rolling resistance (of a single wheel) 
 F(1)= 0.01*Fz*1000 - r(3);         % r(3) = Fr
 Fd= 0.5*Cd*ro*Af*v^2;              %calculating the drag force  
 
 if(sat==0)
    %Fy=M*(r*cos(phi))/4;
    %Fx=Fxmax*(1-(Fy/Fymax)^2)^0.5 ;    %calculating the longitudinal tire force based on elliptical function assumption (see https://www.researchgate.net/figure/The-friction-ellipse-showing-maximum-lateral-and-longitudinal-forces-the-resultant_fig1_224381490)                                      
    F(2)= (Fxmax^2)*(1 -  (((M*r(1)*cos(phi))/4)/Fymax)^2) - (r(2))^2 ; % r(2) = Fx   r(1) = r
 else   
    % when (Fx-Fr)>0.5*TmaxM/Rw
    F(2)= 0.25*TmaxB/Rw + r(3) -r(2) ;  %Brake saturated condition
 end
 %ax=(2*Fx-4*Fr-Fd)/M;
 F(3) = 4*r(2) - (4*r(3) + Fd ) - r(1)*M*sin(phi) ;

