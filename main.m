clear
clc
s=tf('s');
%Some test Functions 
G1= (7*s+1) / ( 7*s^3+21*s^2-28*s );
G2=(s+4)/( s*(s+2)*(s+7)^2 );
G3=(4*s^2 + 23*s + 15)/ ( (s-1)^2*(s^2 + 9) );
tic
nyquistfull(G3,1e-3,1000);
%This function is faster than nyquist()
%blue is for positive frequencies
%red is for negative frequencies
%green is for frequencies around imaginary poles
toc