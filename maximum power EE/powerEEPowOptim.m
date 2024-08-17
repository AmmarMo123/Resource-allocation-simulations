clc; clear all; close all;
addpath('../model');

%Variable initialization
W = 100 * 10^6;              % Total Bandwidth
RCS = 0.1;                   % Target radar cross section
f = 10^9;                    % Carrier frequency
c = 3 * 10^8;                % Speed of light
wavelength = (3 * 10^8) / f; % Wavelength
PathLossRadar = 4.5;         % Path loss exponent for radar
PathLossComm = 2.5;          % Path loss exponent for communication
Pmax = 39.8;                 % Maximum power 30dMb is 1 watt
k_B = 1.380649 * 10^(-23);   % Boltzman Constant
Ttemp = 724;                 % Absolute temperature
gamma1 = 1/3;                % Performance priority
gamma2 = 1/3;                % Performance priority
gamma3 = 1/3;                % Performance priority
m = 3;                       % Nakagami coefficient
J = 1;                       % Number of clutters
g4 = 0.01;                   % Small scale fading of clutter 1
g5 = 0.001;                  % Small scale fading of clutter 2
G = 1;                       % Antenna gain
Rc = (20*10^6) / W ;         % Communication QoS requirement
Rr = (5*10^6) / W;           % Sensing QoS requirement
trials = 3000;
x_values = [];
y_values = [];
unwanted_x = [];
endIndex = 0;
unwanted = 0;
totalUnwanted = 0;
tau = 1 / 3;
different = 0;
circuitPower = 2;
delta = 0.001;
counter = 1;

% Define the center and radius of the circle
radius = 100;
center = [0; 0];

for Pmax = 20:5:45
	total = 0;
	unwanted = 0;
	for i = 1:trials
    	% Creating model, retreiving distances, and seperating the user and clutter distances
    	allDist = createModel(radius, center, J);
    	userDist = allDist(1:3);
    	clutterDist = allDist(4: end);

    	% Preallocate the vector
    	gj = zeros(J, 1);
    	gj(1) = g4;
    	% Fill the vector with values
    	for k = 2:J
        	gj(k) = g5;
    	end
   	 
    	g1 = gamrnd(m, 1/m) * gamrnd(m, 1/m);
    	g2 = gamrnd(m, 1/m) * gamrnd(m, 1/m);
    	h3 = gamrnd(m, 1/m);
       	 
    	a1 = G * userDist(1)^(-PathLossRadar) * RCS * wavelength^2 * g1;
    	a2 = 4 * pi * G * userDist(2)^(-PathLossComm) * c^2 * g2;
    	a3 = G * userDist(3)^(-PathLossComm) * c^2 * h3;
    	b1 = sum(G * clutterDist.^(-PathLossRadar) * RCS * wavelength^2 .* gj);
    	b2 = sum(f^2 * G * clutterDist.^(-PathLossRadar) * RCS * wavelength^2 .* gj);
    	b3 = tau * (4 * pi * f)^2 * k_B * Ttemp * W;
    	c1 = tau * (4 * pi)^3 * k_B * Ttemp * W;
    	c2 = tau * f^2 * (4 * pi)^3 * k_B * Ttemp * W;

    	A1 = a1 / b1;
    	B1 = c1 / b1;
    	A2 = a2 / b2;
    	B2 = c2 / b2;
    	A3 = a3 / b3;

    	counter = 1;
    	F = [100];
    	lambda = [0];
   	 
    	while F(counter)>=delta && counter <= 20
     	cvx_begin
     	cvx_solver Mosek
         	variable y(3) %power
         	maximize((tau * log(1+A1-A1*B1*inv_pos(B1+y(1))))/ log(2) + (tau * log(1+A2-A2*B2*inv_pos(B2+y(2))))/ log(2) + (tau * log(1+ (A3 * y(3)))) / log(2) - lambda(counter) * (circuitPower + y(1) + y(2) + y(3)))
         	subject to
            	(tau * log(1+A1-A1*B1*inv_pos(B1+y(1))))/ log(2) >= Rr;
            	(tau * log(1+A2-A2*B2*inv_pos(B2+y(2))))/ log(2) >= Rc + Rr;
            	(tau * log(1+ (A3 * y(3)))) / log(2) >= Rc;
            	y(1) + y(2) + y(3) <= Pmax;
            	y >= 0;
     	cvx_end
     	if isinf(cvx_optval) || isnan(cvx_optval)
        	break;
     	end
     	counter=counter+1
     	F(counter)=cvx_optval
     	lambda(counter)=((tau * log(1+A1-A1*B1*inv_pos(B1+y(1))))/ log(2) + (tau * log(1+A2-A2*B2*inv_pos(B2+y(2))))/ log(2) + (tau * log(1+ (A3 * y(3)))) / log(2))/(circuitPower + y(1) + y(2) + y(3))
    	end %end of while loop

     	if isinf(cvx_optval) || isnan(cvx_optval)
        	unwanted = unwanted + 1;
        	totalUnwanted = totalUnwanted + 1;
     	else
        	total = total + lambda(counter);
        end
	end
	endIndex = endIndex + 1;
	x_values(endIndex) = Pmax;
	y_values(endIndex) = total / (trials - unwanted);
end

plot(x_values , y_values * W, '-o');
xlabel('Maximum power (watt)');  % Set the x-axis label
ylabel('Energy Efficiency');  % Set the y-axis label