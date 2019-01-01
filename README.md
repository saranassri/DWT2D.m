# DWT2D.m
2D-DWT for the Geophysical data fusion

%% Synthesic data 

%% 2D-DWT based Guo and Yan (2017) algorithm for fused gravity & magnetic data

%%Code by: Sara Nasri and Amin Roshandel kahoo

%%Tel.:+989124395972 E-mail adress:roshandel@shahroodut.ac.ir

%Load Magnetic and Gravity data (Magnetic is image_new_1 and Gravity data is image_new_2)

load image_new_1

load image_new_2

%%Normalize

F=F./max(abs(F(:))); 

E=E./max(abs(E(:)));

%%selection the desired section

x1(:,:)=F(:,35,:); 

x2(:,:)=E(:,35,:);

%%transpose matrix data

x1=x1'; 

x2=x2';

%% made meshgrid for interpolation data

[X,Y]=meshgrid(1:70,1:30); 

[XX,YY]=meshgrid(1:.05:70,1:.01:30);

%%interpolation final models

image_new_1=interp2(X,Y,x1,XX,YY); 

image_new_2=interp2(X,Y,x2,XX,YY);

%%Normalize

image_new_1=image_new_1./max(abs(image_new_1(:))); 

image_new_2=image_new_2./max(abs(image_new_2(:)));

%%Convert matrix to grayscale image

image_new_1=mat2gray(image_new_1); 

image_new_2=mat2gray(image_new_2);

%2-D wavelet decompose of input data in 3-level and got off it approximate and details coefficients (DWT)

[c1,s1]=wavedec2(image_new_1,3,'db1'); 

[c2,s2]=wavedec2(image_new_2,3,'db1');

A1 = appcoef2(c1,s1,'db1',3); 

A2 = appcoef2(c2,s2,'db1',3);

A=.5*A1+.5*A2; 

for i=1:3

    [H1,V1,D1]=detcoef2('all',c1,s1,i);
    
    [H2,V2,D2]=detcoef2('all',c2,s2,i);
    
    H{i}=fusion_detail_2D(H1,H2,1);
    
    V{i}=fusion_detail_2D(V1,V2,1);
    
    D{i}=fusion_detail_2D(D1,D2,1);
    
end

cfusion=[A(:)']; 

for i=3:-1:1

    temp1=H{i};
    
    temp2=V{i};
    
    temp3=D{i};
    
    cfusion=[cfusion temp1(:)' temp2(:)' temp3(:)'];
    
end

%Multilevel 2-D wavelet reconstruction the using inversion of the discrete wavelet transform (IDWT)

xfusion = waverec2(cfusion,s1,'db1'); 
