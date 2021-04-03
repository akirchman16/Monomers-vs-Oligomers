clear all;
close all;

% This simulation will be an ODE model of monomers in solution forming
% tetramers, trimers, dimers and going back until equilibrium. This is not 
% a stochastic model, but rather a deterministic model that uses an Euler 
% method to numerically solve the ODEs. This will be advanced to a model
% that has maximum of general k-mers.

% kinetic rate constants
k_l = 1;    %lengthening of oligomers (units: 1/Ms)
k_s = 1;    %shortening of oligomers (units: 1/s) 

% Initial concentrations of each oligomer (units: M)
L_1(1) = 10;    %monomers
L_2(1) = 0;    %dimers
L_3(1) = 0;    %trimers
L_4(1) = 0;    %tetramers

t(1) = 0;   %initial time

% Values needed to start the system
Loops = 0;
Equilibrium_1 = 0;
Equilibrium_2 = 0;
Equilibrium_3 = 0;
Equilibrium_4 = 0;

% Time step - VARIABLE (units: s)
dt = 0.001;

% Actual system that does all the work
while ~Equilibrium_1 && ~Equilibrium_2 && ~Equilibrium_3 && ~Equilibrium_4
%     Tracks the number of iterations of the code
    Loops = Loops+1;
    i = Loops;
    
%     Euler Method to numerically solve coupled ODEs
    Matrix = [(1-k_l*(L_1(i)+L_2(i)+L_3(i))*dt) (k_s*dt) (k_s*dt) (k_s*dt); (k_l*(L_1(i)-L_2(i))*dt) (1-(k_l*L_2(i)+k_s)*dt) (k_s*dt) (k_s*dt); (k_l*(L_2(i)-L_3(i))*dt) (0) (1-k_s*dt) (k_s*dt); (k_l*L_3(i)*dt) (k_l*L_2(i)*dt) (0) (1-2*k_s*dt)]; %matrix of 'changes'
    State_1 = [L_1(i); L_2(i); L_3(i); L_4(i)]; %matrix of current states
    State_2 = Matrix*State_1; %matrix multiplication to calculate new concentrations
    
%     Updating Concentrations and Time
    L_1(i+1) = State_2(1,1);
    L_2(i+1) = State_2(2,1);
    L_3(i+1) = State_2(3,1);
    L_4(i+1) = State_2(4,1);
    t(i+1) = t(i)+dt;
    
%     Testing for Equilibrium;
    if Loops > 100
        Avg_Change1 = mean(abs(diff(L_1(i-98:i+1))));
        Limits1 = [L_1(i-98),L_1(i+1)];
        Avg_Change2 = mean(abs(diff(L_2(i-98:i+1))));
        Limits2 = [L_2(i-98),L_2(i+1)];
        Avg_Change3 = mean(abs(diff(L_3(i-98:i+1))));
        Limits3 = [L_3(i-98),L_3(i+1)];
        Avg_Change4 = mean(abs(diff(L_4(i-98:i+1))));
        Limits4 = [L_4(i-98),L_4(i+1)];
        if (Avg_Change1 <= 1e-5) && (abs(diff(Limits1)) <= 1e-5)
            Equilibrium_1 = 1;
        end
        if (Avg_Change2 <= 1e-5) && (abs(diff(Limits2)) <= 1e-5)
            Equilibrium_2 = 1;
        end
        if (Avg_Change3 <= 1e-5) && (abs(diff(Limits3)) <= 1e-5)
            Equilibrium_3 = 1;
        end
        if (Avg_Change4 <= 1e-5) && (abs(diff(Limits4)) <= 1e-5)
            Equilibrium_4 = 1;
        end
    end
end

figure(1);
scatter(t,L_1,3,'r','filled');
hold on;
scatter(t,L_2,3,'b','filled');
scatter(t,L_3,3,'g','filled');
scatter(t,L_4,3,'k','filled');
xlabel('Time (s)');
xlim([0 max(t)]);
ylabel('Concentrations (M)');
ylim([min([min(L_1),min(L_2),min(L_3),min(L_4)]) max([max(L_1),max(L_2),max(L_3),max(L_4)])]);
title('Monomers, Dimers, and Trimers in Solution');
legend('Monomers','Dimers','Trimers','Tetramers');