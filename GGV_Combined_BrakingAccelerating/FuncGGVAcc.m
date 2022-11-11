function F = FuncGGVAcc(r,v,phi,sat)
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
    F(2)= 0.5*TmaxM/Rw + r(3) -r(2) ;   %Motor saturated condition
 end
 %ax=(2*Fx-4*Fr-Fd)/M;
 F(3) = 2*r(2) - (4*r(3) + Fd ) - r(1)*M*sin(phi) ;

