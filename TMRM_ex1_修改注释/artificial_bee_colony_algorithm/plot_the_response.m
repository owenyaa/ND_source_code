%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:32
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 1 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}%
%          	\label{eq4.2}%
%          	\begin{cases}%
%          		\begin{aligned}%
%          		  & \dot{S}(t)+S^2(t)=1, t\geq 0 \\%
%          		  & S(0)=0                       %
%          		\end{aligned}%
%          	\end{cases}%
%          \end{equation}. The structure of the program is as follows:%
%          The basis function in this program is exponential function \chi_i(t)=e^{-it}%
%          The MATLAB codes about the artificial bee_colony (ABC) algorithm%
%% To plot the response of the ABC results and the exact solutions

clear;clc;close all;
order=49;Tdata=0:0.01:200;
load 'matlab.mat';
a0=GlabolParams(1,1);
GlabolParams1=GlabolParams';
disx=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
disx(1,:)=a0;
for i=1:order
    disx(1,:)=disx(1,:)+GlabolParams1(i+1,1).*exp(-i.*Tdata);
    dx(1,:)=dx(1,:)-i*GlabolParams1(i+1,1).*exp(-i.*Tdata);
end
y=tanh(Tdata);
%figure;
plot(Tdata,disx(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,y,'r-','MarkerSize',8);
h1=legend('$$ABC$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
