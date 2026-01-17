% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775
% ΚΕΡΑΙΕΣ ΚΑΙ ΔΙΑΔΟΣΗ ΕΡΓΑΣΙΑ 2 - ΘΕΜΑ 2

clc; clear; close all;

%% ΟΡΙΣΜΟΣ ΠΑΡΑΜΕΤΡΩΝ ΣΧΕΔΙΑΣΗΣ
f0 = 400;                   % Κεντρική συχνότητα λειτουργίας (MHz)
lambda = 300 / f0;          % Μήκος κύματος (m) -> 0.75m για 400MHz
N_turns = 10;               % Αριθμός σπειρών

% Διαστάσεις Έλικας (Axial Mode)
C = lambda;                 % Περιφέρεια ίση με το μήκος κύματος
R_helix = C / (2*pi);       % Ακτίνα της έλικας
S = lambda / 4;             % Βήμα έλικας (Spacing) -> λ/4
Wire_Radius = lambda/200;   % Ακτίνα σύρματος (Διάμετρος ~ λ/100)

% Παράμετροι Επιπέδου Γείωσης (Ground Plane)
R_ground = lambda / 2;      % Συνολική ακτίνα ground plane (λ/2)
Gap = lambda / 40;          % Απόσταση τροφοδοσίας από το ground
r_rings = [0.25, 0.50, 0.75, 1.00] * R_ground; % Ακτίνες 4 ομόκεντρων δακτυλίων

%% ΔΗΜΙΟΥΡΓΙΑ ΑΡΧΕΙΟΥ NEC
filename = 'helix_geometry.nec';
fid = fopen(filename, 'w');
fprintf(fid, 'CM Helix Antenna with Spider-web Ground Plane\n');
fprintf(fid, 'CM Designed by Karatis Dimitrios (10775) - f0=400MHz\n');
fprintf(fid, 'CE\n');

%% ΚΑΤΑΣΚΕΥΗ GROUND PLANE (RADIALS & RINGS)
% --- Α. Ειδική Ακτίνα στις 0 μοίρες (Σημείο Σύνδεσης Τροφοδοσίας) ---
nodes_0 = unique([0, R_helix, r_rings]);
nodes_0 = sort(nodes_0);
for i = 1:length(nodes_0)-1
    start_r = nodes_0(i);
    end_r = nodes_0(i+1);
    fprintf(fid, 'GW  1  1  %.4f  0  0  %.4f  0  0  %.4f\n', ...
            start_r, end_r, Wire_Radius);
end

% --- Β. Υπόλοιπες 7 Ακτίνες (Ανά 45 μοίρες) ---
nodes_others = [0, r_rings]; 
for k = 1:7
    angle = k * 45;
    for m = 1:length(nodes_others)-1
        r_start = nodes_others(m);
        r_end = nodes_others(m+1);
        x1 = r_start * cosd(angle); y1 = r_start * sind(angle);
        x2 = r_end * cosd(angle);   y2 = r_end * sind(angle);
        fprintf(fid, 'GW  1  1  %.4f  %.4f  0  %.4f  %.4f  0  %.4f\n', ...
                x1, y1, x2, y2, Wire_Radius);
    end
end

% --- Γ. Ομόκεντροι Δακτύλιοι (Spider-web Rings) ---
for r = r_rings
    for k = 0:7
        ang1 = k * 45; ang2 = (k+1) * 45;
        x1 = r * cosd(ang1); y1 = r * sind(ang1);
        x2 = r * cosd(ang2); y2 = r * sind(ang2);
        fprintf(fid, 'GW  2  3  %.4f  %.4f  0  %.4f  %.4f  0  %.4f\n', ...
                x1, y1, x2, y2, Wire_Radius);
    end
end

%% ΣΥΡΜΑ ΤΡΟΦΟΔΟΣΙΑΣ (FEED WIRE)
fprintf(fid, 'GW  100  1  %.4f  0  0  %.4f  0  %.4f  %.4f\n', ...
        R_helix, R_helix, Gap, Wire_Radius);

%% ΚΑΤΑΣΚΕΥΗ ΕΛΙΚΑΣ (SEGMENT-BY-SEGMENT)
Segs_per_turn = 17;
Total_segs = N_turns * Segs_per_turn;
d_theta = 2*pi / Segs_per_turn;
d_z = S / Segs_per_turn;
x_prev = R_helix; y_prev = 0; z_prev = Gap;
for i = 1:Total_segs
    theta = i * d_theta;
    z_curr = Gap + i * d_z;
    x_curr = R_helix * cos(theta);
    y_curr = R_helix * sin(theta);
    fprintf(fid, 'GW  4  1  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
            x_prev, y_prev, z_prev, x_curr, y_curr, z_curr, Wire_Radius);
    x_prev = x_curr; y_prev = y_curr; z_prev = z_curr;
end

fprintf(fid, 'GE  0\n');

%% ΠΗΓΗ ΔΙΕΓΕΡΣΗΣ ΚΑΙ ΣΥΧΝΟΤΗΤΑ
% Τοποθέτηση πηγής τάσης (Voltage Source) στο σύρμα τροφοδοσίας (Tag 100)
fprintf(fid, 'EX  0  100  1  0  1.0  0.0\n');

% Ορισμός Κεντρικής Συχνότητας (FR) - Μία συχνότητα χωρίς sweep
% FR 0 (Τύπος) | 1 (Αριθμός συχνοτήτων) | 0 | 0 | f0 (MHz) | 0 (Step)
fprintf(fid, 'FR  0  1  0  0  %.2f  0\n', f0);

% Τέλος αρχείου
fprintf(fid, 'EN\n');
fclose(fid);

disp('Το αρχείο helix_geometry.nec δημιουργήθηκε επιτυχώς!');