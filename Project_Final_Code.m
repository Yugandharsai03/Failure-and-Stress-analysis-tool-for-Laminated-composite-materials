disp('Considering the Analysis for Graphite Polymer Composite');
%The given properties are for Graphite Polymer Composite.
E1=155;E2=12.10;E3=12.10;%E Values
V23=0.458;V13=0.248;V12=0.248;V21=0.01936;%V Values
G23=3.20;G13=4.40;G12=4.40;%G Values
alpha1=-0.01800*10^-6;alpha2=24.3*10^-6;alpha3=24.3*10^-6;%alpha values
beta1=146*10^-6;beta2=4770*10^-6;beta3=4770*10^-6;%beta values
%Defining the Angle Properties.
layers=input('Enter the number of layers: ');
a=zeros(1,layers);
for value=1:layers
    a(value)=input('Enter the angle of orientation for the Ply: ');
value=value+1;
end
a;%Angle
value=zeros(1,layers);
for value=1:layers
    value(value)=cosd(a(value));
    value=value+1;
end
value
m=zeros(1,layers);
for value=1:layers
    m(value)=cosd(a(value));
    value=value+1;
end
m

n=zeros(1,layers);
for value=1:layers
    n(value)=sind(a(value));
    value=value+1;
end
n
%Thickness
t=150*10^-6;
T=layers*t
disp('T')

%Obtaining the Compliance Matrix.
S11=1/E1,S12=-V12/E1,S13=-V13/E1,
S21=-V12/E1,S22=1/E2,S23=-V23/E2,
S31=-V13/E1,S32=-V23/E2,S33=1/E3,
S44=1/G23,S55=1/G13,S66=1/G12;
S=[S11,S12,S13,0,0,0;S21,S22,S23,0,0,0;S31,S32,S33,0,0,0;0,0,0,S44,0,0;0,0,0,0,S55,0;0,0,0,0,0,S66]%Compliance Matrix
disp('Compliance Matrix is')
S

%Obtaining the Stiffness Matrix.
S_=(S11*S22*S33)-(S11*S23*S23)-(S22*S13*S13)-(S33*S12*S12)+2*(S12*S23*S13)
C11=((S22*S33)-(S23*S23))/S_,C12=((S13*S23)-(S12*S33))/S_,C22=((S33*S11)-(S13*S13))/S_,C13=((S12*S23)-(S13*S22))/S_,
C33=((S11*S22)-(S12*S12))/S_,C23=((S12*S13)-(S23*S11))/S_,C44=1/S44,C55=1/S55,C66=1/S66
C=[C11,C12,C13,0,0,0;C12,C22,C23,0,0,0;C13,C23,C33,0,0,0;0,0,0,C44,0,0;0,0,0,0,C55,0;0,0,0,0,0,C66]%Stiffness Matrix
disp('Stiffness Matrix is')
C

%Obtaining Reduced Compliance Matrix.
R_C=[S11,S12,0;S12,S22,0;0,0,S66]%Reduced Compliance Matrix.
disp('Reduced Compliance Matrix is')
R_C

%Obtaining Reduced Stiffness Matrix.
Q11=C11-(C13)^2/C33,Q12=C12-(C13*C23)/C33,Q22=C22-(C23)^2/C33,Q66=C66
R_K=[Q11,Q12,0;Q12,Q22,0;0,0,Q66]%Reduced Stiffness Matrix.
disp('Reduced Stiffness Matrix is')
R_K
%Obtaining Transformed Reduced Compliance.
for value=1:layers
    S_11=S11*m(value)^4+((2*S12)+(S66))*n(value)^2*m(value)^2+S22*n(value)^4
    value=value+1;
end
for value=1:layers
    S_12=(S11+S22-S66)*n(value)^2*m(value)^2+(S12)*n(value)^4+m(value)^4
    value=value+1;
end
for value=1:layers
    S_16=(2*S11-2*S12-S66)*n(value)*m(value)^3-(2*S22-2*S12-S66)*n(value)^3*m(value)
    value=value+1;
end
for value=1:layers
    S_22=(S11)*n(value)^4+(2*S12+S66)*n(value)^2*m(value)^2+S22*m(value)^4
    value=value+1;
end
for value=1:layers
    S_26=(2*S11-2*S12-S66)*n(value)^3*m(value)-(2*S22-2*S12-S66)*n(value)*m(value)^3
    value=value+1;
end
for value=1:layers
    S_66=2*(2*S11+2*S22-4*S12-S66)*n(value)^2*m(value)^2+S66*((n(value)^4)+(m(value)^4))
    value=value+1;
end
T_R_C=[S_11,S_12,S_16;S_12,S_22,S_26;S_16,S_26,S_66]%Transformed Reduced Compliance
disp('Transformed Reduced Compliance Matrix is')
t_r_c=transpose(T_R_C)%Calculating Transposed Matrix

%Obtaining Transformed Reduced Stiffness.
for value=1:layers
    Q_11(value)=Q11*m(value)^4+((2*Q12)+2*(Q66))*n(value)^2*m(value)^2+Q22*n(value)^4
    value=value+1;
end
for value=1:layers
    Q_12(value)=(Q11+Q22-4*Q66)*n(value)^2*m(value)^2+(Q12)*n(value)^4+m(value)^4
    value=value+1;
end

for value=1:layers
    Q_16(value)=(Q11-Q12-2*Q66)*n(value)*m(value)^3+(Q12-Q22+2*Q66)*n(value)^3*m(value)
    value=value+1;
end
for value=1:layers
    Q_22(value)=(Q11)*n(value)^4+(2*Q12+2*Q66)*n(value)^2*m(value)^2+Q22*m(value)^4
    value=value+1;
end
for value=1:layers
    Q_26(value)=(Q11-Q12-2*Q66)*n(value)^3*m(value)+(Q12-Q22+2*Q66)*n(value)*m(value)^3
    value=value+1;
end
for value=1:layers
    Q_66(value)=(Q11+Q22-2*Q12-2*Q66)*n(value)^2*m(value)^2+Q66*((n(value)^4)+(m(value)^4))
    value=value+1;
end
T_R_S=[Q_11,Q_12,Q_16;Q_12,Q_22,Q_26;Q_16,Q_26,Q_66] %Transformed Reduced Stiffness
t_r_S=[sum(Q_11);sum(Q_12);sum(Q_16);sum(Q_22);sum(Q_26);sum(Q_66)];
t_r_s=transpose(t_r_S);%Calculating Transposed Matrix
disp('Transformed Reduced Stiffness is: ')
T_R_S

%Obtaining Engineering Properties in Global Coordinate System. 
disp('Engineering Properties of the each Lamina in the Global Coordinate system: ')
 
for value= 1:layers
    E_x(value)=(E1)/((m(value)^4+((E1/G12)-2*V12)*n(value)^2*m(value)^2)+((E1/E2))*n(value)^4);
    value=value+1;    
end
E_x
 
 for value= 1:layers
    V_xy(value)=(V12*(n(value)^4 + m(value)^4)-((1+(E1/E2)-(E1/G12)*n(value)^2*m(value)^2)))/(m(value)^4 + ((E1/G12)-2*V12)*n(value)^2*m(value)^2 + (E1/E2)*n(value)^4);
    value=value+1;    
 end  
 V_xy
for value=1:layers
    E_y(value)=(E2)/((m(value)^4+(((E2/G12)-2*V12)*n(value)^2*m(value)^2)+(E2/E1)*n(value)^4));
    value=value+1;
end
E_y

for value=1:layers
    V_yx(value)=((V21*(n(value)^4+m(value)^4))-(1+(E2/E1)-(E2/G12)*n(value)^2*m(value)^2))/((m(value)^4)+(((E2/G12)-2*V21)*n(value)^2*m(value)^2)+((E2/E1)*n(value)^4));
    value = value+1;
end
V_yx
for value=1:layers
    G_xy(value)=(G12)/(n(value)^4+m(value)^4+2*(2*(G12/E1)*(1+2*V12)+2*(G12/E2)-1)*n(value)^2*m(value)^2); 
    value = value+1;
end
G_xy
%To Find the A Matrix.

A = zeros(1,layers);
for value = 1:6
    A(value) = sum(t_r_S(value))*T;
   value= value+1;
end
 disp('A-Matrix')
A
 
 
%To Find the B Matrix.
    
B = zeros(1,layers);
for value = 1:6
    B(value) = (1/2)*(sum(t_r_S*(value))*(T*(value+1)^2-T*(value)^2));
   value= value+1;
end
 disp('B-Matrix')
B
 
%To Find the D Matrix.

D = zeros(1,layers);
for value = 1:6
    D(value) = (1/3)*(sum(t_r_S*(value))*(T*(value+1)^3-T*(value)^3));
   value= value+1;
end
 disp('D-matrix')
D

%To Find the ABD matrix.
 
 disp('The ABD matrix')
ABD =[A(1) A(2) A(3) B(1) B(2) B(3);
    A(2) A(4) A(5) B(2) B(4) B(5);
    A(3) A(5) A(6) B(3) B(5) B(6);
    B(1) B(2) B(3) D(1) D(2) D(3);
    B(2) B(4) B(5) D(2) D(4) D(5);
    B(3) B(5) B(6) D(3) D(5) D(6)]
  
% abd Matrix
 disp('The abd matrix')
 abd = inv(ABD)

%Effective Engineering Properties for Laminate
fprintf('Effective Engineering Properties for Laminate: ')
    

for value=1:layers
    G_xy = A(6)./T;
    value=value+1;
end
G_xy
for value=1:layers
   V_xy = A(2)./A(4);
    value=value+1;
end
V_xy
for value=1:layers
  V_yx = A(2)./A(1);
    value=value+1;
end
V_yx
for value=1:layers
 E_x = (A(1)*A(4)- (A(2)^2))./(A(4)*T);
    value=value+1;
end
 E_x
for value=1:layers
 E_y = (A(1)*A(4)- (A(2)^2))./(A(1)*T);
    value=value+1;
end
E_y




 
 %By using the abd matrix we can obtain Laminate's extensional strains and
 %curvatures of the laminate.
 
 % To find strains:   
strain_xo = input('extensional strain x_direction ');strain_yo = input('extensional strain y_direction ');strain_xy0 = input('extensional strain x_y direction ');curvature_xo = input('curvature x_ direction ');curvature_yo= input('curvature y_ direction ');curvature_xyo= input('curvature x_y direction ');
strain_x = zeros(1,layers);strain_y = zeros(1,layers);strain_xy = zeros(1,layers);

for value = 1:layers    
strain_x(value) = ((strain_xo)  + (a(value))*(curvature_xo));
strain_y(value) = ((strain_yo)  + (a(value))*(curvature_yo));
strain_xy(value)= ((strain_xy0) + (a(value))*(curvature_xyo));
value = value+1;
end
strain = [strain_x; strain_y; strain_xy]
 
%To find stresses:
for value=1:layers
     stress_1(value) = [Q_11(value).*strain_x(value) + Q_12(value).*strain_y(value) + Q_16(value).*strain_xy(value)]; 
     value=value+1;
end
 
for value= 1:layers
    stress_2(value) = [Q_12(value).*strain_x(value) + Q_22(value).*strain_y(value) + Q_26(value).*strain_xy(value)]; 
    value = value+1;
end
for value= 1:layers
    stress_3(value) = [Q_16(value).*strain_x(value) + Q_26(value).*strain_y(value) + Q_66(value).*strain_xy(value)];
    value = value+1;
end
 Stress = [stress_1; stress_2; stress_3]
 
 %T-Matrix Representation:
 %T_matrix= [(value).^2,n.^2,2*(value)*n;
            %n.^2,(value).^2,-2*(value)*n;
           %(-value)*n,(value)*n,(value).^2-n.^2;]
        
%To find the principle stresses The T-Matrix is multiplied by axial stresses:

for value= 1:layers
    principle_stress_1(value) = [m(value)^2.*stress_1(value)+ n(value)^2.*stress_2(value) + 2*m(value)*n(value)*stress_3(value)];
    value = value+1;
end
 
for value= 1:layers
    principle_stress_2(value) = [n(value)^2.*stress_1(value)+ m(value)^2.*stress_2(value) - 2*m(value)*n(value)*stress_3(value)];
    value = value+1;
end
 
for value= 1:layers
    principle_stress_3(value) = [-m(value).*n(value).*stress_1(value) + m(value).*n(value).*stress_2(value) + (m(value)^2-n(value)^2).*stress_3(value)];
    value = value+1;
end
 
principle_stress = [A; B; C]
 
% Maximum Stress Criteria.
sigma1c = (-1.250),sigma1t = (1.50),sigma2c = (-0.200),sigma2t = (0.050),tow12 = (0.100);
A=principle_stress_1;B=sigma1c;C=sigma1t;D=principle_stress_2;E=sigma2c;F=sigma2t;
G=principle_stress_3;H=-tow12;I=tow12;
%Finding analysis for Tension failure criteria.
for value = 1:layers
if abs(A(value)) >= abs(B) || abs(A(value)) <= abs(C)
   fprintf('Safe___Tension \n',value);
else
   fprintf('Failure___Tension \n',value);
end
value = value+1;
end
%Finding analysis for Compression failure criteria.  
for value = 1:layers
if abs(D(value)) >= abs(E) || abs(D(value)) <= abs(F)
    fprintf('Safe___Compression \n',value);
else
    fprintf('Failure___Compression \n',value);
end
value = value+1;
end
 %Finding analysis for Shear Failure Criteria. 
for value = 1:layers
if abs(G(value)) >= abs(H) || abs(G(value)) <= abs(I)
    fprintf('Safe___Shear \n',value);
else
    fprintf('Failure___Shear \n',value);
end
value = value+1;
end
 
%Analysis of T-SaiWu Failure Criteria. 
for value=1:layers
 
    if abs(((1/C)+(1/B))*principle_stress_1(value)) + (((1/F)+(1/E))*principle_stress_2(value)) + ((-1/(C*B))*principle_stress_1(value)^2) + ((-1/(F*E))*principle_stress_2(value)^2) + ((1/(I))^2)*principle_stress_3(value)^2 + (2*((-0.5)*sqrt((-1/(C*B))*(-1/(F*E))*principle_stress_1(value)*principle_stress_2(value)))) <= 1
           fprintf('Safe___T-Sai Wu \n',value);
    else
           fprintf('Failure___T-Sai Wu \n',value);
    end
end
%Analysis of T-Sai Hill Criteria.
for value=1:layers
    if abs(((principle_stress_1(value))/(1/C)+(1/B))^2+(((principle_stress_2(value))/(1/F)+(1/E))^2+(I)/(0))^2-((principle_stress_1(value))*((principle_stress_2(value))/((1/C)+(1/B))^2))) < 1
        fprintf('Safe___T-Sai Hill \n',value);
    else
        fprintf('Failure___T-Sai Hill \n',value);
    end
end