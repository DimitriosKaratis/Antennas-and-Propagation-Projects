% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775
% ΚΕΡΑΙΕΣ ΚΑΙ ΔΙΑΔΟΣΗ ΕΡΓΑΣΙΑ 2 - ΘΕΜΑ 1

clc; clear; close all;

% Παράμετροι βάσει λ = 0.6m 
lambda = 0.6; 
Rd = 0.34 * lambda;            % Ακτίνα δίσκου 
Lc = 0.5 * lambda;             % Μήκος συρμάτων κώνου 
theta0 = deg2rad(30);          % Μισή γωνία ανοίγματος (2theta0=60) 
wire_rad = (lambda/200) / 2;   % Ακτίνα συρμάτων 
gap = lambda/20;               % Απόσταση τροφοδοσίας 
n_wires = 8;                   % Αριθμός συρμάτων  
n_seg = 19;                    % Κατάτμηση 

% Δημιουργία αρχείου .nec
fid = fopen('discone_geometry.nec', 'w');
fprintf(fid, 'CM Discone Antenna Geometry Template\n');
fprintf(fid, 'CM Designed for lambda = 0.6m\n');
fprintf(fid, 'CE\n');

tag = 1;

% --- ΓΕΩΜΕΤΡΙΑ ΔΙΣΚΟΥ ---
for i = 1:n_wires
    phi = (i-1) * (2*pi/n_wires);
    x2 = Rd * cos(phi);
    y2 = Rd * sin(phi);
    % GW [tag] [segments] [x1] [y1] [z1] [x2] [y2] [z2] [radius] 
    fprintf(fid, 'GW %d %d 0 0 0 %f %f 0 %f\n', tag, n_seg, x2, y2, wire_rad);
    tag = tag + 1;
end

% --- ΓΕΩΜΕΤΡΙΑ ΚΩΝΟΥ ---
z_cone_top = -gap;
for i = 1:n_wires
    phi = (i-1) * (2*pi/n_wires);
    x_b = Lc * sin(theta0) * cos(phi);
    y_b = Lc * sin(theta0) * sin(phi);
    z_b = z_cone_top - (Lc * cos(theta0));
    fprintf(fid, 'GW %d %d 0 0 %f %f %f %f %f\n', tag, n_seg, z_cone_top, x_b, y_b, z_b, wire_rad);
    tag = tag + 1;
end

% --- ΣΥΡΜΑ ΤΡΟΦΟΔΟΣΙΑΣ (FEED WIRE) ---
% Ενώνει το κέντρο του δίσκου (0,0,0) με την κορυφή του κώνου (0,0,-gap) 
fprintf(fid, 'GW %d 1 0 0 0 0 0 %f %f\n', tag, z_cone_top, wire_rad);

fprintf(fid, 'GE 0\n'); % Τέλος Γεωμετρίας

fprintf(fid, 'EX 0 %d 1 0 1 0\n', tag); % Τροφοδοσία στο τελευταίο σύρμα 
fprintf(fid, 'FR 0 1 0 0 500\n');       % Default συχνότητα 500 MHz (f0 = 300/0.6) 
fprintf(fid, 'EN\n');

fclose(fid);
disp('Το αρχείο discone_geometry.nec δημιουργήθηκε επιτυχώς!');