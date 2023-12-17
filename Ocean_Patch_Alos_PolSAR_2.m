clc;
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D:\Himanshu_Laptop\Himanshu_Laptop_Drive\PolSAR Data\IEEE GRSL\Urban Patch Mumbai\U3_SIG\T3\T11 --- enter the file location address where you store T11 element

% fid1 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T11.bin");
% fid2 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T22.bin");
% fid3 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T22.bin");
% fid4 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T12_real.bin");
% fid5 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T12_imag.bin");
% fid6 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T13_real.bin");
% fid7 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T13_imag.bin");
% fid8 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T23_real.bin");
% fid9 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\SIG\T3\T23_imag.bin");


%%%%%%%%%%%%% Urban5 patch %%%%%%%%%%%%%%%%%%%

fid1 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T11.bin");
fid2 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T22.bin");
fid3 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T33.bin");
fid4 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T12_real.bin");
fid5 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T12_imag.bin");
fid6 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T13_real.bin");
fid7 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T13_imag.bin");
fid8 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T23_real.bin");
fid9 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_urban5\SIG\T3\T23_imag.bin");

Nligfin=290;
Ncol=260;

%%%%%%%%%%%%% Ocean Patch %%%%%%%%%%%%%%%%%%
% fid1 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T11.bin");
% fid2 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T22.bin");
% fid3 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T33.bin");
% fid4 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T12_real.bin");
% fid5 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T12_imag.bin");
% fid6 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T13_real.bin");
% fid7 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T13_imag.bin");
% fid8 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T23_real.bin");
% fid9 = fopen("C:\Users\panka\Downloads\Alos-2 San Francisco Patches\Alos_San_ocean1\T3\T23_imag.bin");
% 
% Nligfin=217;
% Ncol=229;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for lig = 1 : Nligfin
    T11(lig,1:Ncol) = fread(fid1, Ncol, 'float32');
    T22(lig,1:Ncol) = fread(fid2, Ncol, 'float32');
    T33(lig,1:Ncol) = fread(fid3, Ncol, 'float32');
    T12_real(lig,1:Ncol) = fread(fid4, Ncol, 'float32');
    T12_imag(lig,1:Ncol) = fread(fid5, Ncol, 'float32');
    T13_real(lig,1:Ncol) = fread(fid6, Ncol, 'float32');
    T13_imag(lig,1:Ncol) = fread(fid7, Ncol, 'float32');
    T23_real(lig,1:Ncol) = fread(fid8, Ncol, 'float32');
    T23_imag(lig,1:Ncol) = fread(fid9, Ncol, 'float32'); 
    %C13_real(lig,1:Ncol) = fread(fid10, Ncol, 'float32'); 
end


c=T11-T22;
span=T11+T22+T33;
T12=complex(T12_real,T12_imag);
T13=complex(T13_real,T13_imag);
T23=complex(T23_real,T23_imag);
T21=complex(T12_real, -T12_imag);
T31=complex(T13_real, -T13_imag);
T32=complex(T23_real,-T23_imag);

T_v=1/4*[2 0 0; 0 1 0; 0 0 1];
for i=1:Nligfin
    for j=1:Ncol
        %T=[T11(i,j) complex(T12_real(i,j),T12_imag(i,j)) complex(T13_real(i,j),T13_imag(i,j));...
          %conj(complex(T12_real(i,j),T12_imag(i,j))) T22(i,j) complex(T23_real(i,j),T23_imag(i,j));...
          %conj(complex(T13_real(i,j),T13_imag(i,j))) conj(complex(T23_real(i,j),T23_imag(i,j))) T33(i,j)];
        T=[T11(i,j) T12(i,j) T13(i,j); T21(i,j) T22(i,j)  T23(i,j);  T31(i,j) T32(i,j)  T33(i,j)];
        T(isnan(T))=0;
        eig_val = eig(T,T_v); %it will have 3 eigenvalues Lambda_1,Lambda_2, Lambda_3.
        pv1(i,j) = min(eig_val);  %Here the volume scattering power is equal to the minimum of the eigenvalue
    end
end

T11r=T11-0.5*pv1;
T22r=T22-0.25*pv1;
T33r=T33-0.25*pv1;
% for i=1:Nligfin
%     for j=1:Ncol
%         Tr=[T11r(i,j) T12(i,j) T13(i,j); T21(i,j) T22r(i,j) T23(i,j); T31(i,j) T32(i,j) T33r(i,j)];
%         thetafp1(i,j)=atand((T11r(i,j)+T22r(i,j)+T33r(i,j))*(T11r(i,j)-T22r(i,j)-T33r(i,j))/(T11r(i,j)*(T22r(i,j)+T33r(i,j))+ (T11r(i,j)+T22r(i,j)+T33r(i,j))*(T11r(i,j)+T22r(i,j)+T33r(i,j))));
%         ps1(i,j)=(T11r(i,j)+T22r(i,j)+T33r(i,j))*(1+sin(2*thetafp1(i,j)))/2;
%         pd1(i,j)=(T11r(i,j)+T22r(i,j)+T33r(i,j))*(1-sin(2*thetafp1(i,j)))/2;
%     end
% end
tic; 
Tr=[T11r T12 T13; T21 T22r T23; T31 T32 T33r];

%%%%%%%%%%%%%%%%%%%%%%%% OAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1=(0.25)*atand(2*T12_real./(T11r-T22r));
%thetafp1=atand((T11r+T22r+T33r)*(T11r-T22r-T33r)/(T11r*(T22r+T33r)+ (T11r+T22r+T33r)*(T11r+T22r+T33r)));
%ps1=(T11r+T22r+T33r)*(1+sind(2*thetafp1))/2;
%pd1=(T11r+T22r+T33r)*(1-sind(2*thetafp1))/2;
T11n=T11r.*(cos(2*v1).*cos(2*v1))+ T22r.*(sin(2*v1).*sin(2*v1))+(T12_real.*sin(4*v1));
T22n=T22r.*(cos(2*v1).*cos(2*v1))+ T11r.*(sin(2*v1).*sin(2*v1))-(T12_real.*sin(4*v1));
T12n=T12_real.*(cos(4*v1)) + 1i*(T12_imag) +(T22r-T11r).*sin(4*v1)/2;
T13n=T13.*cos(2*v1) + T23.*(sin(2*v1));
T23n=T23.*cos(2*v1)- T13.*(sin(2*v1));
T33n=T33r;

v2=(0.25)*atan(2*T12_imag./(T11n-T22n));
T11n1=T11n.*(cos(2*v2).*cos(2*v2))+ T22.*(sin(2*v2).*sin(2*v2))+(T12_imag.*sin(4*v2));
T22n1=T22n.*(cos(2*v2).*cos(2*v2))+ T11.*(sin(2*v2).*sin(2*v2))-(T12_imag.*sin(4*v2));
T12n1=T12_imag.*(cos(4*v2)) + 1i*(T22-T11).*sin(4*v2)/2;
T13n1=T13.*cos(2*v2) + 1i*T23.*(sin(2*v2));
T23n1=T23.*cos(2*v2) + 1i*T13.*(sin(2*v2));
T33n1=T33;
T32n1=conj(T23n1);
T31n1=conj(T13n1);
% T_n1=[T11n1 T12n1 T13n1; T21n1 T22n1 T23n1; T31n1 T32n1 T33n1];
% deter_Tn1=det(T_n1);
% trace_Tn1=trace(T_n1);
% span_Tn1=trace(T_n1);
% m_fp_n1=sqrt(1-((27*deter_Tn1)/trace_Tn1^3));

for i=1:Nligfin
    for j=1:Ncol   
        T_n1=[T11n1(i,j) 0 T13n1(i,j); 0 T22n1(i,j) T23n1(i,j) ;  T31n1(i,j) T32n1(i,j)  T33n1(i,j)];
%         deter_Tn1=det(T_n1);
        trace_Tn1=trace(T_n1);
%         span_Tn1=trace(T_n1);
%         m_fp_n1=sqrt(1-((27*deter_Tn1)/trace_Tn1^3));
        thetafp2(i,j)=real(atand(trace_Tn1*(T11n1(i,j)-T22n1(i,j)-T33n1(i,j))/(T11n1(i,j)*(T22n1(i,j)+T33n1(i,j))+(trace_Tn1*trace_Tn1))));
        ps1(i,j)=real(trace_Tn1)*(1+sind(2*thetafp2(i,j))/2);
        pd1(i,j)=real(trace_Tn1)*(1-sind(2*thetafp2(i,j))/2);
    end
end
% 
% x = real((0.25)*atand(2*real(T12_real)./(T11r - T22r)));
% for i=1:Nligfin
%     for j=1:Ncol
%         TR=[T11r(i,j) complex(T12_real(i,j),T12_imag(i,j)) complex(T13_real(i,j),T13_imag(i,j));...
%           conj(complex(T12_real(i,j),T12_imag(i,j))) T22r(i,j) complex(T23_real(i,j),T23_imag(i,j));...
%           conj(complex(T13_real(i,j),T13_imag(i,j))) conj(complex(T23_real(i,j),T23_imag(i,j))) T33r(i,j)];
%         R_x=[cosd(2*x(i,j)) sind(2*x(i,j)) 0; -sind(2*x(i,j)) cosd(2*x(i,j)) 0; 0 0 1];
%         T_x=R_x*TR*R_x';
%         y = real((0.25)*atand(2*imag(T_x(1,2))./(T_x(1,1) - T_x(2,2))));
%         R_y = [cosd(2*y) 1i*sind(2*y) 0; 1i*sind(2*y) cosd(2*y) 0; 0 0 1];
%         T_y = R_y*T_x*R_y';        
%         deter_T_y = det(T_y);
%         tra_T_y = trace(T_y);
%         span_T_y = trace(T_y);
%         m_FP_T_y = sqrt(1 - ((27*deter_T_y)/tra_T_y^3));
%         theta_FP_T_y = real(atand((m_FP_T_y*span_T_y*(T_y(1,1) - T_y(2,2) - T_y(3,3)))/((T_y(1,1)*(T_y(2,2) + T_y(3,3)))+((m_FP_T_y^2)*(span_T_y^2)))));   
%         ps1(i,j) = real(trace(T_y)*(1 + sind(2*theta_FP_T_y))/2);
%         pd1(i,j) = real(trace(T_y)*(1 - sind(2*theta_FP_T_y))/2);
%     end
% end

for i=1:Nligfin
    for j=1:Ncol
        if ps1(i,j)>pd1(i,j) && ps1(i,j)>pv1(i,j)       %% zone 8 blue%%
           psm1(i,j,1)=0;
           psm1(i,j,2)=0;
           psm1(i,j,3)=1;
        elseif pd1(i,j)>=ps1(i,j) && pd1(i,j)>=pv1(i,j)       %% zone 8 red%%
           psm1(i,j,1)=1;
           psm1(i,j,2)=0;
           psm1(i,j,3)=0;
        elseif pv1(i,j)>ps1(i,j) && pv1(i,j)>pd1(i,j)   %% zone 6 green%%
           psm1(i,j,1)=0;
           psm1(i,j,2)=1;
           psm1(i,j,3)=0;

        end
    end
end

s_ps1=0;
s_pd1=0;
s_span=0;
s_pv1=0;


% 
ps1=ps1./span;
pd1=(pd1)./span;
pv1=pv1./span;

for i=1:Nligfin
    for j=1:Ncol
        s_ps1=ps1(i,j)+s_ps1;
        s_pd1=pd1(i,j)+s_pd1;
        s_pv1=pv1(i,j)+s_pv1;
    end
end

m_ps1=s_ps1/(Nligfin*Ncol);
m_pd1=s_pd1/(Nligfin*Ncol);
m_pv1=s_pv1/(Nligfin*Ncol);
m_span=m_ps1+m_pd1+m_pv1;




