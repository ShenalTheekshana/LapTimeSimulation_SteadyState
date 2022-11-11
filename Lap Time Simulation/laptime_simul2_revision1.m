%Loading the GGV Diagram Data
load GGV5_Acc_Brake.mat

figure(1)
[c,h] = contour(X,Y,Z,290);

contourTable = getContourLineCoordinates(c);

%% Importing curvature values
c = xlsread('Optimum LapTrack.xlsx','Sheet3','B2:B971'); %importing from the excel flie
s = xlsread('Optimum LapTrack.xlsx','Sheet3','C2:C971'); %importing from the excel flie

%% Extracting values from GGV Diagram
%{
v_current = v;
[~,idx]=min(abs(contourTable.Level - v_current));
v_closest = contourTable.Level(idx);
matrix_Range = find(contourTable.Level==v_closest);
x_Matrix = contourTable.X(matrix_Range);
y_Matrix = contourTable.Y(matrix_Range);

x_Matrix_acc = [];
y_Matrix_acc = [];
x_Matrix_brk = [];
y_Matrix_brk = [];

[max_x_val,max_x_index] = max(x_Matrix);

for k = 1:length(matrix_Range)
    if y_Matrix(k) < y_Matrix(max_x_index)
        x_Matrix_brk = [x_Matrix_brk;x_Matrix(k)];
        y_Matrix_brk = [y_Matrix_brk;y_Matrix(k)];
    else
        x_Matrix_acc = [x_Matrix_acc;x_Matrix(k)];
        y_Matrix_acc = [y_Matrix_acc;y_Matrix(k)];
    end
end

a_lat_max = max_x_val;

[max_y_acc_val,max_y_acc_index] = max(y_Matrix_acc);
a_lon_max = max_y_acc_val;

[max_y_brk_val,max_y_brk_index] = min(y_Matrix_brk);
brk_lon_max = abs(max_y_brk_val);
%}

%% Initiating Matrices  

%Matirix to store exiting velocity of eac segment while accelerating
v_out_acc = [];
%Matirix to store exiting velocity of eac segment while breaking
v_out_brk = [];
lon_acc_mat_a = [];
lon_acc_mat_b = [];

lat_acc_mat_a = [];
lat_acc_mat_b = [];

%% 1ST SEGMENT: Initial velociy is zero.

v_in = 0;
%interpolating GGV Diagram
a_lat_max     =   ;
a_lon_max     =   ;
brk_lon_max     =   ;
%--------------------------------------------------------------------------------
v_current = v_in;
[val,idx]=min(abs(contourTable.Level - v_current));
v_closest = contourTable.Level(idx);
matrix_Range = find(contourTable.Level==v_closest);
x_Matrix = contourTable.X(matrix_Range);
y_Matrix = contourTable.Y(matrix_Range);

x_Matrix_acc = [];
y_Matrix_acc = [];
x_Matrix_brk = [];
y_Matrix_brk = [];

[max_x_val,max_x_index] = max(x_Matrix);

for k = 1:length(matrix_Range)
    if y_Matrix(k) < y_Matrix(max_x_index)
        x_Matrix_brk = [x_Matrix_brk;x_Matrix(k)];
        y_Matrix_brk = [y_Matrix_brk;y_Matrix(k)];
    else
        x_Matrix_acc = [x_Matrix_acc;x_Matrix(k)];
        y_Matrix_acc = [y_Matrix_acc;y_Matrix(k)];
    end
end

a_lat_max = max_x_val;

[max_y_acc_val,max_y_acc_index] = max(y_Matrix_acc);
a_lon_max = max_y_acc_val;

%---------------------------------------------------------------------------------

v_max = sqrt(abs(a_lat_max / c(1))); %c is the matrix containing curvatures of each segment
v_out = sqrt(v_in^2 + 2 * a_lon_max * s(1)); %s is the matix containing arc lengths of each segment

if v_out > v_max
    v_out = v_max
        
    lon_acc = 0;
    lat_acc = a_lat_max;
else
    lon_acc = a_lon;
    lat_acc = a_lat;
end

v_out_acc = [v_out_acc;v_out];
lon_acc_mat_a = [lon_acc_mat_a;lon_acc];
lat_acc_mat_a = [lat_acc_mat_a;lat_acc];

%% Acceleration (Any other segment)

for i=2:length(s)
    
    v_in = v_out;
    % corresponding to 'v_in' value by GGV diagram
    %a_lat_max     =   ;
    %a_lon_max     =   ;
%-----------------------------------------------------------------------------------------------------------  
    v_current = v_in;
    [val,idx]=min(abs(contourTable.Level - v_current));
    v_closest = contourTable.Level(idx);
    matrix_Range = find(contourTable.Level==v_closest);
    x_Matrix = contourTable.X(matrix_Range);
    y_Matrix = contourTable.Y(matrix_Range);

    x_Matrix_acc = [];
    y_Matrix_acc = [];
    x_Matrix_brk = [];
    y_Matrix_brk = [];

    [max_x_val,max_x_index] = max(x_Matrix);

    for k = 1:length(matrix_Range)
        if y_Matrix(k) < y_Matrix(max_x_index) %dividing the matrix into two, longitudinal acc and logitudinal brk.
            x_Matrix_brk = [x_Matrix_brk;x_Matrix(k)];
            y_Matrix_brk = [y_Matrix_brk;y_Matrix(k)];
        else
            x_Matrix_acc = [x_Matrix_acc;x_Matrix(k)];
            y_Matrix_acc = [y_Matrix_acc;y_Matrix(k)];
        end
    end

    a_lat_max = max_x_val;

    [max_y_acc_val,max_y_acc_index] = max(y_Matrix_acc);
    a_lon_max = max_y_acc_val;

    [max_y_brk_val,max_y_brk_index] = min(y_Matrix_brk);
    brk_lon_max = abs(max_y_brk_val);
%------------------------------------------------------------------------------------------------------------  
    
  
    if c(i) ~= 0
        v_max = sqrt(abs(a_lat_max / c(i)));
    else
        v_max = 30;
    end
    
    if v_max < v_in
        v_in = v_max
        v_out = v_in
    else
        a_lat = v_in^2 * c(i)
        %a_lon = sin(acos(a_lat / a_lat_max)).* a_lat_max; %this is considering GGV diagram is cylindrical(should be changed)
        %Possible a_lon while maintaining a_lat should be found
        a_lon = abs(interp1(x_Matrix_acc,y_Matrix_acc,a_lat,'spline'));
        
        v_out = sqrt(v_in^2 + 2 * a_lon * s(i));
    end
    
    if v_out > v_max
        v_out = v_max
        
        lon_acc = 0;
        lat_acc = a_lat_max;
    else
        lon_acc = a_lon;
        lat_acc = a_lat;
    end
    
    v_out_acc = [v_out_acc;v_out];
    lon_acc_mat_a = [lon_acc_mat_a;lon_acc];
    lat_acc_mat_a = [lat_acc_mat_a;lat_acc];
        
    
end

%% Breaking (Any other segment) 

v_in = sqrt(abs(a_lat_max / c(length(c))));

v_out_brk = [v_out_brk;v_in]; 

for i=flip(2:length(s))
    
    v_out = v_in;
    % corresponding to 'v_in' value by GGV diagram
    %a_lat_max     =   ; 
    %brk_lon_max    =   ;
%---------------------------------------------------------------------------------------------------------------    
    v_current = v_out;
    [val,idx]=min(abs(contourTable.Level - v_current));
    v_closest = contourTable.Level(idx);
    matrix_Range = find(contourTable.Level==v_closest);
    x_Matrix = contourTable.X(matrix_Range);
    y_Matrix = contourTable.Y(matrix_Range);

    x_Matrix_acc = [];
    y_Matrix_acc = [];
    x_Matrix_brk = [];
    y_Matrix_brk = [];

    [max_x_val,max_x_index] = max(x_Matrix);

    for k = 1:length(matrix_Range)
        if y_Matrix(k) < y_Matrix(max_x_index)
            x_Matrix_brk = [x_Matrix_brk;x_Matrix(k)];
            y_Matrix_brk = [y_Matrix_brk;y_Matrix(k)];
        else
            x_Matrix_acc = [x_Matrix_acc;x_Matrix(k)];
            y_Matrix_acc = [y_Matrix_acc;y_Matrix(k)];
        end
    end

    a_lat_max = max_x_val;

    [max_y_acc_val,max_y_acc_index] = max(y_Matrix_acc);
    a_lon_max = max_y_acc_val;

    [max_y_brk_val,max_y_brk_index] = min(y_Matrix_brk);
    brk_lon_max = abs(max_y_brk_val);
%-----------------------------------------------------------------------------------------------------------
    
    if c(i) ~= 0
        v_max = sqrt(abs(a_lat_max / c(i)));
    else
        v_max = 30;
    end    
    
    if v_max < v_out
        v_out = v_max
        v_in = v_out
    else
        a_lat = v_out^2 * c(i)
        %Possible b_lon while maintaining the  a_lat should be found
        brk_lon = abs(interp1(x_Matrix_brk,y_Matrix_brk,a_lat,'spline')); 
        
        v_in = sqrt(v_out^2 + 2 * brk_lon * s(i));
    end
    
    if v_in > v_max
        v_in = v_max;
        
        lon_acc = 0;
        lat_acc = a_lat_max;
    else
        lon_acc = brk_lon;
        lat_acc = a_lat;
    end
    
    v_out_brk = [v_out_brk;v_in];
    lon_acc_mat_b = [lon_acc_mat_b;lon_acc];
    lat_acc_mat_b = [lat_acc_mat_b;lat_acc];
    
        
    
end
v_out_brk = flip(v_out_brk)
lon_acc_mat_b = flip(lon_acc_mat_b)
lat_acc_mat_b = flip(lat_acc_mat_b)

%% Final Velocity

v_final = min(v_out_acc,v_out_brk);

%% Exporting Data

writematrix(v_out_acc,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','E2:E971');
writematrix(v_out_brk,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','F2:F971');
writematrix(v_final,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','G2:G971');
writematrix(lat_acc_mat_a,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','I2:I971');
writematrix(lat_acc_mat_b,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','J2:J971');
writematrix(lon_acc_mat_a,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','L2:L971');
writematrix(lon_acc_mat_b,'Optimum LapTrack.xlsx','Sheet','Sheet3','Range','M2:M971');

%{
figure(1)
plot(d,v_out_acc,'-r');
title('Acceration')
figure(2)
plot(d,v_out_brk,'-b');
title('Deceleration')
figure(3)
plot(d,v_final,'-g');
title('Final')
figure(4)
plot(d,v_out_acc,'-r');
hold on
plot(d,v_out_brk,'-b');
hold on
plot(d,v_final,'-g');
hold on
title('All-together')
figure(5)
plot(d,c(1:928,1))
title('Curvature')
figure(3)
plot(d,v_final,'-g','LineWidth',2);
title('Final')
 %}
    
