%Loading the GGV Diagram Data
load GGV5_Acc_Brake.mat

figure(1)
[c,h] = contour(X,Y,Z,290);

contourTable = getContourLineCoordinates(c);

%% Importing curvature values
c = xlsread('LapDistanceCurvatureData.xlsx','Sheet1','B2:B929'); %importing from the excel flie
s = xlsread('LapDistanceCurvatureData.xlsx','Sheet2','C2:C929'); %importing from the excel flie

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

%% 1ST SEGMENT: Initial velociy is zero.

v_in = 0;
%interpolating GGV Diagram
%a_lat_max     =   ;
%a_lon_max     =   ;
%brk_lon_max     =   ;
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
end

v_out_acc = [v_out_acc;v_out]; 

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
        
    end
    
    v_out_acc = [v_out_acc;v_out]; 
        
    
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
        v_in = v_max
        
    end
    
    v_out_brk = [v_out_brk;v_in]; 
        
    
end
v_out_brk = flip(v_out_brk)

%% Final Velocity

v_final = min(v_out_acc,v_out_brk);
v_final_kmph = v_final.*(3600/1000);



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


