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
trials = 1000;
x_values = [];
y_values = [];
radarMI = [];
isacMI = [];
commMI = [];
unwanted_values = [];
endIndex = 0;
unwanted = 0;
totalUnwanted = 0;
Gamma3 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.66];

% Define the center and radius of the circle
radius = 100;
center = [0; 0];

for j = 1:length(Gamma3)
    total = 0;
    Gamma1 = (1 - Gamma3(j))/2;
    Gamma2 = (1 - Gamma3(j))/2;
    unwanted = 0;
    radarMITotal = 0;
    isacMITotal = 0;
    commMITotal = 0;
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
        
        g1 = randraw('nakagami', [m, 1/m], 1) *  randraw('nakagami', [m, 1/m], 1);
        g2 = randraw('nakagami', [m, 1/m], 1) *  randraw('nakagami', [m, 1/m], 1);
        g3 = randraw('nakagami', [m, 1/m], 1) *  randraw('nakagami', [m, 1/m], 1);

        a1 = 4 * pi * G * userDist(1)^(-PathLossComm) * c^2 * g1;
        a2 = 4 * pi * G * userDist(2)^(-PathLossComm) * c^2 * g2;
        a3 = 4 * pi * G * userDist(3)^(-PathLossComm) * c^2 * g3;
        b1 = sum(f^2 * G * clutterDist.^(-PathLossRadar) * RCS * wavelength^2 .* gj);
        b2 = sum(f^2 * G * clutterDist.^(-PathLossRadar) * RCS * wavelength^2 .* gj);
        b3 = sum(f^2 * G * clutterDist.^(-PathLossRadar) * RCS * wavelength^2 .* gj);
        c1 = f^2 * (4 * pi)^3 * k_B * Ttemp * W;
        c2 = f^2 * (4 * pi)^3 * k_B * Ttemp * W;
        c3 = f^2 * (4 * pi)^3 * k_B * Ttemp * W;
        
        A1 = (b1/c1) + (a1/c1);
        C1 = (b1/c1);
        A2 = (b2/c2) + (a2/c2);
        C2 = (b2/c2);
        A3 = (b3/c3) + (a3/c3);
        C3 = (b3/c3);

        cvx_begin
        cvx_solver Mosek
            variable y(3) %P
            variable x(3) %tau
            maximize((Gamma1 * (1/(1-A1/C1)*(rel_entr(A1*y(1)+x(1),C1*y(1)+x(1))+A1/C1*rel_entr(C1*y(1)+x(1),A1*y(1)+x(1)))))/log(2) + (Gamma2 * (1/(1-A2/C2)*(rel_entr(A2*y(2)+x(2),C2*y(2)+x(2))+A2/C2*rel_entr(C2*y(2)+x(2),A2*y(2)+x(2))))) / log(2) + (Gamma3(j) * (1/(1-A3/C3)*(rel_entr(A3*y(3)+x(3),C3*y(3)+x(3))+A3/C3*rel_entr(C3*y(3)+x(3),A3*y(3)+x(3))))) / log(2))
            subject to
               (1/(1-A1/C1)*(rel_entr(A1*y(1)+x(1),C1*y(1)+x(1))+A1/C1*rel_entr(C1*y(1)+x(1),A1*y(1)+x(1)))) / log(2) >= Rr;
               (1/(1-A2/C2)*(rel_entr(A2*y(2)+x(2),C2*y(2)+x(2))+A2/C2*rel_entr(C2*y(2)+x(2),A2*y(2)+x(2)))) / log(2) >= Rc + Rr;
               (1/(1-A3/C3)*(rel_entr(A3*y(3)+x(3),C3*y(3)+x(3))+A3/C3*rel_entr(C3*y(3)+x(3),A3*y(3)+x(3)))) / log(2) >= Rc;
                y(1) + y(2) + y(3) <= Pmax;
                x(1) + x(2) + x(3) == 1;
                y >= 0;
                x >= 0; 
        cvx_end

        if isinf(cvx_optval) || isnan(cvx_optval)
           unwanted = unwanted + 1;
           totalUnwanted = totalUnwanted + 1;
        else
           total = total + cvx_optval;
           radarMITotal = radarMITotal + (x(1) * log2(1+(a1 * y(1))/((b1*y(1)) + (c1*x(1)))));
           isacMITotal = isacMITotal + (x(2) * log2(1+(a2 * y(2))/((b2*y(2)) + (c2*x(2)))));
           commMITotal = commMITotal + (x(3) * log2(1+(a3 * y(3))/((b3*y(3)) + (c3*x(3)))));
        end
    end
    endIndex = endIndex + 1;
    x_values(endIndex) = Gamma3(j);
    y_values(endIndex) = total / (trials - unwanted);
    radarMI(endIndex) = radarMITotal / (trials - unwanted);
    isacMI(endIndex) = isacMITotal / (trials - unwanted);
    commMI(endIndex) = commMITotal / (trials - unwanted);
    unwanted_values(endIndex) = unwanted;
end

plot(x_values , radarMI, '-o', 'LineWidth', 3);
hold on
plot(x_values , isacMI, '-o', 'LineWidth', 3);
plot(x_values , commMI, '-o', 'LineWidth', 3);
legend('Radar MI', 'ISaC MI', 'Comm Rate');
xlabel('Priority of communication only user, Gamma_3');  % Set the x-axis label
ylabel('MI or data rate (bps)');  % Set the y-axis label
hold off