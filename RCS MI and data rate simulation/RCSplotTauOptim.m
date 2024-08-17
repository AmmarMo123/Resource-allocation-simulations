clc; clear all; close all;
addpath('../model');

%Variable initialization
W = 100 * 10^6;              % Total Bandwidth
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
RCS = [0.1, 0.3, 0.5, 1, 2, 5, 10]; % Target radar cross section
trials = 3000;
x_values = [];
y_values = [];
endIndex = 0;
unwanted = 0;
totalUnwanted = 0;
total0 = 0;
P = Pmax / 3;

% Define the center and radius of the circle
radius = 100;
center = [0; 0];

for j = 1:length(RCS)
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
            
        a1 = P * G * userDist(1)^(-PathLossRadar) * RCS(j) * wavelength^2 * g1;
        a2 = P * 4 * pi * G * userDist(2)^(-PathLossComm) * c^2 * g2;
        a3 = P * G * userDist(3)^(-PathLossComm) * c^2 * h3;
        b1 = P * sum(G * clutterDist.^(-PathLossRadar) * RCS(j) * wavelength^2 .* gj);
        b2 = P * sum(f^2 * G * clutterDist.^(-PathLossRadar) * RCS(j) * wavelength^2 .* gj);
        b3 = (4 * pi * f)^2 * k_B * Ttemp * W;
        c1 = (4 * pi)^3 * k_B * Ttemp * W;
        c2 = f^2 * (4 * pi)^3 * k_B * Ttemp * W;
        
        y1 = a1 / c1;
        A1 = b1 / c1;
        y2 = a2 / c2;
        A2 = b2 / c2;
        A3 = a3 / b3;
         
        cvx_begin
        cvx_solver Mosek
            variable x(3) %tau
            maximize((-rel_entr(x(1)+A1,x(1)+A1+y1) - A1*(rel_entr((x(1)+A1)/y1,(x(1)+A1)/y1+1) + rel_entr((x(1)+A1)/y1+1,(x(1)+A1)/y1)))/log(2) + (-rel_entr(x(2)+A2,x(2)+A2+y2) - A2*(rel_entr((x(2)+A2)/y2,(x(2)+A2)/y2+1) + rel_entr((x(2)+A2)/y2+1,(x(2)+A2)/y2)))/log(2) + -rel_entr(x(3),x(3)+A3)/ log(2)) 
            subject to
                (-rel_entr(x(1)+A1,x(1)+A1+y1) - A1*(rel_entr((x(1)+A1)/y1,(x(1)+A1)/y1+1) + rel_entr((x(1)+A1)/y1+1,(x(1)+A1)/y1)))/log(2)  >= Rr;
                (-rel_entr(x(2)+A2,x(2)+A2+y2) - A2*(rel_entr((x(2)+A2)/y2,(x(2)+A2)/y2+1) + rel_entr((x(2)+A2)/y2+1,(x(2)+A2)/y2)))/log(2) >= Rc + Rr;
                -rel_entr(x(3),x(3)+A3) / log(2) >= Rc;
                x(1) + x(2) + x(3) == 1;
                x >= 1e-6; 
        cvx_end

        if isinf(cvx_optval) || isnan(cvx_optval)
           unwanted = unwanted + 1;
           totalUnwanted = totalUnwanted + 1;
        else
           total = total + cvx_optval;
        end
    end
    endIndex = endIndex + 1;
    x_values(endIndex) = RCS(j);
    y_values(endIndex) = total / (trials-unwanted);
end

plot(x_values , y_values * W, '-o');
xlabel('Target radar cross section (m^2)');  % Set the x-axis label
ylabel('Sum of MI and data rate (bps)');  % Set the y-axis label