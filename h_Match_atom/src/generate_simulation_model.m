function atom_sim = generate_simulation_model()
%% cell distance, bond length and layer distance of primitive unitcell
latticeC = 2.7;                                                            % lattice constant unite (Angstrom)
bondL = latticeC/sqrt(3);                                                  % bondlength unite (Angstrom)
layerD = 3.4;                                                              % layer distance unite (Angstrom)

%% three vector for the carbon graphene Hexagonal
a1 = [1*latticeC, 0*bondL];
a2 = [latticeC/2, sqrt(3)*latticeC/2];

%% the matrix used to shift the tpye 1 atom to type 2 atom
un_pos=[0, 0; 0, bondL];

%% generate the basic Hexagonal graphene 
Row = -45:45;                                                              % the row number
Col = -25:80;                                                              % the column number
Base = zeros(length(Row),length(Col),2);                                   % (row, col, x-y position)
for tempi = 1:length(Row)
    for tempj = 1:length(Col)
        Base(tempi,tempj,:) = Row(tempi).*a1 + Col(tempj).*a2;
    end
end
Base_2D = reshape(Base,length(Row)*length(Col),2);
clear tempi tempj a1 a2 bondL latticeC

%% generate the first layer
Pos = [un_pos(1,:)+Base_2D;-un_pos(2,:)+Base_2D];
clear Base Base_2D un_pos Col Row

%% generate the bilayer graphene
% angle m  , n          ;  angle  m , n       
% 0.100 331, 330 1310764;  1.696  20, 19 4564
% 0.200 166, 165 328684 ;  2.005  17, 16 3268 
% 0.300 111, 110 146524 ;  2.134  16, 15 2884
% 0.400 83 , 82  81676  ;  2.281  15, 14 2524
% 0.497 67 , 66  53068  ;  2.450  14, 13 2188
% 0.596 56 , 55  36964  ;  2.646  13, 12 1876
% 0.797 42 , 41  20668  ;  2.876  12, 11 1588
% 0.987 34 , 33  13468  ;  3.150  11, 10 1324
% 1.018 33 , 32  12676  ;  3.481  10, 9  1084
% 1.050 32 , 31  11908  ;  3.890  9 , 8  868
% 1.085 31 , 30  11164  ;  4.408  8 , 7  676
% 1.121 30 , 29  10444  ;  5.086  7 , 6  508
% 1.539 22 , 21  5548   ;  6.009  6 , 5  364
% 1.614 21 , 20  5044   ;  7.341  5 , 4  244
% 1.696 20 , 19  4564   ;  9.430  4 , 3  148
% 1.788 19 , 18  4108   ;  13.174 3 , 2  76
% 1.890 18 , 17  3676   ;  21.787 2 , 1  28

%% calculate the twist angle
m = 17;
n = 16;
twist_angle=acosd(1/2*(m^2+n^2+4*m*n)/(m^2+n^2+m*n));
% rotation matrix
% anticlockwise rotate
R1 = [cosd(twist_angle/2) -sind(twist_angle/2); sind(twist_angle/2) cosd(twist_angle/2)];
% clockwise rotate
R2 = [cosd(-twist_angle/2) -sind(-twist_angle/2); sind(-twist_angle/2) cosd(-twist_angle/2)];

% anticlockwise rotate
Layer1=(R1*Pos')';Layer1(:,3) = +layerD./2;
% anticlockwise rotate
Layer2=(R2*Pos')';Layer2(:,3) = -layerD./2;
atom_sim = [Layer1',Layer2'];
clear Layer1 Layer2 layerD Pos twist_angle
end