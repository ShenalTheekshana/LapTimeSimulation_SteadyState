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

opts = optimoptions('fsolve', 'TolFun', 1E-3, 'TolX', 1E-3,'MaxFunctionEvaluations',1000000,'MaxIterations',1000000);
for i = 0:3:30
    v = i;      %defining velocity (z in cylindrical coordinates)
    %Braking (FuncGGVBraking used)
    for j = -16:1:0
        phi = j/10;   % defining phi (azimuth angle in cylindrical coordinates)
        x = fsolve(@(x)FuncGGVBraking(x,v,phi,0),[1,0,0], opts); %solving the function "FuncGGVBraking" by inputing... 
                                                                 %v,phi and variable to switch from normal condition... 
                                                                 %to brake saturated (maximum brake torque reached)condition
        R(i/3+1,j+17) = x(1);       %writing r (radius of cylindrical coorinates) to matrix R
        
        if((x(2)+x(3))>0.25*TmaxB/Rw)                                   %Checking the if maximum brake torque has reached
             x = fsolve(@(x)FuncGGVBraking(x,v,phi,1),[1,0,0], opts);   %the 4th variable here(sat) is set to 1 to switch... 
                                                                        %to brake saturated condition
             R(i/3+1,j+17) = x(1);  %writing r (radius of cylindrical coorinates) to matrix R
        end
    end
    
    %Accelerating (FuncGGVAcc used)fi
    for j = 0:1:16
        phi = j/10;     % defining phi (azimuth angle in cylindrical coordinates)
        x = fsolve(@(x)FuncGGVAcc(x,v,phi,0),[1,0,0], opts); %solving the function "FuncGGVAcc" by inputing... 
                                                             %v,phi and variable to switch from normal condition... 
                                                              %to motor saturated (maximum motor torque reached)condition
        R(i/3+1,j+17) = x(1);     %writing r (radius of cylindrical coorinates) to matrix R
        
        if((x(2)-x(3))>0.5*TmaxM/Rw)                                    %Checking the if maximum motor torque has reached
             x = fsolve(@(x)FuncGGVAcc(x,v,phi,1),[1,0,0], opts);       %the 4th variable here(sat) is set to 1 to switch... 
                                                                        %to motor saturated condition
             R(i/3+1,j+17) = x(1);   %writing r (radius of cylindrical coorinates) to matrix R
        end
    end

end    
        
v= 0:3:30;                    %creating a vector/array with v values
phi= -1.6:0.1:1.6;            %creating a vector/array with phi values
[PHI,V] = meshgrid(phi,v);    %this creats a mesh grid of v and phi values (google "meshgrid of matlab")
figure(1)
X = R.*cos(PHI);       %converting cylindrical coordinates to cartesian coordinates 
Y = R.*sin(PHI);
Z = V;

figure(2)
surf(X,Y,Z);           %plotting the surface
grid on;
xlabel('ay')
ylabel('ax')
zlabel('v')



