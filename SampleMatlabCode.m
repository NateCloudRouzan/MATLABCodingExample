clear all;%clear all variables
close all;%Close all figures
%Plots are in the bottom section of this code. uncomment/comment a chunk script to
%see related plots.  also save is at the bottom of the script

%Load salmon data from a CSV file. 
load salmon_data.csv;

%Create a time vector that goes the length of the data
t = (1:length(salmon_data)).';

% establish Q to be a symetric matrix that will solve RMS best fit
Q = [sum(t.^2), sum(t);
    sum(t), length(t)];

%Vector of solutions to Q*P = R
R = [sum(t.*salmon_data);
     sum(salmon_data)];

P = Q\R; %solution to Q*P = R

A1 = [Q,R,P]; %Save [Q, R and P into a concatinated matrix]


% %calculate different polonomial coefficients with data 
poly1 = polyfit(t, salmon_data, 1); %first order (Y = Ax + B)
poly2 = polyfit(t, salmon_data, 2); %second order (Y = Ax^2 + bx + C)
poly5 = polyfit(t, salmon_data, 5); %fifth order (Y =Ax^5 + Bx^4 + Cx^3 + Dx^2 + Ex + F)
poly8 = polyfit(t, salmon_data, 8); %eighth order (Y = Ax^8 + Bx^7 + Cx^6 + Dx^5 + Ex^4 + Fx^3 + Gx^2 + Hx + I)

A2 = poly2; %Save second order polynomial in A2 
A3 = poly5; %Save fifth order polynomial in A3
A4 = poly8; %Save eighth orfer polynomial in A4

%Make a prediction for all of these for 2015
year_for_guess = 78; %2015 - 1938

%put the three guesses in A5
A5 = [polyval(poly2, year_for_guess);
      polyval(poly5, year_for_guess);
      polyval(poly8, year_for_guess)];

%Create coarse salmon data so that we can interpolate data
coarse_salmon_data = salmon_data(1:4:end);
t_coarse = t(1:4:end); % save matching time vector

A6 = coarse_salmon_data; %Save the coarse salmon dat to A6 

%take different interpolations of coarse data
 interp_nearest =  interp1(t_coarse, coarse_salmon_data, t,'nearest');
 interp_linear = interp1(t_coarse, coarse_salmon_data, t,'linear');
 interp_cubic = interp1(t_coarse, coarse_salmon_data, t,'cubic');
 interp_spline = interp1(t_coarse, coarse_salmon_data, t,'spline');

 %Save these interpolated data points in a row vectors
 A7 = [interp_nearest, interp_linear, interp_cubic, interp_spline];
 
 %Take RMS errors of different methods and save as A8
 A8 = [sqrt(sum(abs(salmon_data - interp_nearest).^2)/length(t));
       sqrt(sum(abs(salmon_data - interp_linear).^2)/length(t));
       sqrt(sum(abs(salmon_data - interp_cubic).^2)/length(t));
       sqrt(sum(abs(salmon_data - interp_spline).^2)/length(t)); ];

  %time vector defined in lab
  R_t = 0:0.1:1;
  
  %Data given to us in lab
  R =  [-0.088, 4.73, 3.036, -0.299, -1.56, -1.64, -0.466, -0.318, -0.282, -0.289, -0.672]; 
  R1 = [0, 1.25, 2.02, 0.46, -3.74, -7.08, -3.08, 10.80, 24.27, 15.08, -29.70];
  R2 = [0.01, 0.014, 0.022, 0.033, 0.049, 0.073, 0.11, 0.16, 0.24, 0.36, 0.54]; 
  R3 = [0, 10.2, 6.69, 0.63, -2.04, -1.57, -0.27, 0.39, 0.36, 0.09, -0.073];
  R4 = [0.134, 0.146, 0.174, 0.206, 0.236, 0.256, 0.263, 0.256, 0.236, 0.206, 0.174];
  
  %symetric matrix 
  M = [sum(R1.^2), sum(R1.*R2), sum(R1.*R3), sum(R1.*R4)
       sum(R2.*R1), sum(R2.^2), sum(R2.*R3), sum(R2.*R4)
       sum(R3.*R1), sum(R3.*R2), sum(R3.^2), sum(R3.*R4)
       sum(R4.*R1), sum(R4.*R2), sum(R4.*R3), sum(R4.^2)];

  %vector of answers 
  b = [sum(R.*R1)
      sum(R.*R2)
      sum(R.*R3)
      sum(R.*R4)];

  %alpha vector is the solution of M^-1*b
  alpha = M\b;
  
  A9 = M; %Save matrix to A9
  A10 = b; %save vector b to A10
  A11 = alpha; %save vector alpha to A11
  
  %approximate R(t) @ time 10
  A12 = alpha(1)*12.3 + alpha(2)*(-1.45) + alpha(3)*(6.2) + alpha(4)*(-0.03);
  
%Save variables A1-A12  
    save('A1.dat','A1','-ascii');
    save('A2.dat','A2','-ascii');
    save('A3.dat','A3','-ascii');
    save('A4.dat','A4','-ascii');
    save('A5.dat','A5','-ascii');
    save('A6.dat','A6','-ascii');
    save('A7.dat','A7','-ascii');
    save('A8.dat','A8','-ascii');
    save('A9.dat','A9','-ascii');
    save('A10.dat','A10','-ascii');
    save('A11.dat','A11','-ascii');
    save('A12.dat','A12','-ascii');

%figures are below (uncomment block of code to display)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This figure displays the three polynomial approx
%to give the user a better undertanding of the accuracy

plot(t+1937,salmon_data, '.'), hold on;
x_axis = linspace(1938,2014, 1000);
x = linspace(1,77, 1000);
plot(x_axis, polyval(P, x)); %
plot(x_axis,polyval(poly1,x));
plot(x_axis,polyval(poly2,x));
plot(x_axis,polyval(poly5,x));
plot(x_axis,polyval(poly8,x));

%Correct X axis and stuff
    figure(1)
    legend('data points','X^2 ', 'x^5' , 'x^8', 'location', 'best');
    xlabel('Year');
    ylabel('# of Salmon');
    title('Annual Chinook Salmon Count');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This plot displays the interpolations of the course data
%  Uncomment to see this figure
    figure(2)
    plot(t+1937,interp_nearest), hold on;
    plot(t+1937,interp_linear);
    plot(t+1937,interp_cubic);
    plot(t+1937,interp_spline);
    legend('nearest neigbor','linear', 'cubic' , 'spline', 'location', 'best');
    xlabel('Year');
    ylabel('# of Salmon');
    title('Interpolated salmon data');
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This plot displays the four Mori-Zwanzig funct
% in a 2x2plot and also a second figure that displays
% our original function
% (Uncomment to see this plot)
  figure(3);
  subplot(2,2,1);plot(R_t,R1);
    xlabel('time');
    ylabel('R1(t)');
    title('Reduced order R1(t)');
  subplot(2,2,2);plot(R_t,R2);
    xlabel('time');
    ylabel('R2(t)');
    title('Reduced order R2(t)');
  subplot(2,2,3);plot(R_t,R3);
    xlabel('time');
    ylabel('R3(t)');
    title('Reduced order R3(t)');
  subplot(2,2,4);plot(R_t,R4);
    xlabel('time');
    ylabel('R4(t)');
    title('Reduced order R4(t)');
  figure(4);plot(R_t,R) 
       xlabel('time');
        ylabel('R(t)');
        title('base function R(t)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This plot will display how close
%your Mori-Zwanzig is to the function
%(Uncomment to see this plot)
    figure(5)
    plot(R_t,R)
    hold on
    plot(R_t,alpha(1)*R1+alpha(2)*R2+alpha(3)*R3+alpha(4)*R4,'r')
    legend('Actual function', 'approximation', 'location', 'best');
    xlabel('time');
    ylabel('R(t)');
    title('Mori-Zwanzig Approximation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
