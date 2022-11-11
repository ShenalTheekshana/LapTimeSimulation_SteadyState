%Loading the GGV Diagram Data
load GGV5_Acc_Brake.mat

figure(1)
[c,h] = contour(X,Y,Z,290);

contourTable = getContourLineCoordinates(c);
[max_lat_Acc,~] = max(abs(contourTable.X));
[max_lon_Acc,~] = max(contourTable.Y);
[min_lon_Acc,~] = min(contourTable.Y);

Cd_front = 0.003228*(v^3) - 0.001536*(v^2) - 0.008723*v + 0.8459
Cd_Rear = 0.0001991*(v^2) - 0.003648*v + 0.7156
Cd_Under = 0.0001333*(v^3) - 0.009029*(v^2) + 0.1875*v - 0.594

Cl_front = -0.02607*(v^2) + 0.07178*v + 3.641
Cl_Rear = 0.00024*(v^2) + 0.00244*v + 3.016
Cl_Under = -1.098*(v^2) + 0.887*v + 3.411

A_front = 1.7731;
A_rear = 1.6717;
A_under = 1.7336;

Cd = ((Cd_front*A_front)+(Cd_Rear*A_rear)+(Cd_Under*A_under))/(A_front + A_rear + A_under);
Cl = ((Cl_front*A_front)+(Cl_Rear*A_rear)+(Cl_Under*A_under))/(A_front + A_rear + A_under);

Cd = (((0.003228*(v^3) - 0.001536*(v^2) - 0.008723*v + 0.8459)*A_front)+((0.0001991*(v^2) - 0.003648*v + 0.7156)*A_rear)+((0.0001333*(v^3) - 0.009029*(v^2) + 0.1875*v - 0.594
)*A_under))/(A_front + A_rear + A_under);
Cl = ((( -0.02607*(v^2) + 0.07178*v + 3.641)*A_front)+((0.00024*(v^2) + 0.00244*v + 3.016)*A_rear)+((-1.098*(v^2) + 0.887*v + 3.411)*A_under))/(A_front + A_rear + A_under);