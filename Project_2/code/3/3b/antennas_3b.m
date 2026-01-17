% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775
% ΚΕΡΑΙΕΣ ΚΑΙ ΔΙΑΔΟΣΗ ΕΡΓΑΣΙΑ 2 - ΘΕΜΑ 3β

clc; clear; close all;

%% ΟΡΙΣΜΟΣ ΠΑΡΑΜΕΤΡΩΝ
f0 = 12;                    
lambda = 300 / f0;          
L_total = 5 * lambda;       
h = lambda / 4;            

a_wire = (lambda / 200) / 2;   % Κύριος αγωγός
a_feed = a_wire / 4;           % Λεπτά vertical wires (για feed και ground)

segs_hor = 101; 
segs_vert = 5; 

R_load = 317.5;                % Αντίσταση τερματισμού (Ω)

%% ΔΗΜΙΟΥΡΓΙΑ ΑΡΧΕΙΟΥ NEC
filename = 'traveling_wave_geometry.nec';
fid = fopen(filename, 'w');
fprintf(fid, 'CM Traveling Wave Horizontal Antenna (Beverage Type)\n');
fprintf(fid, 'CM Source at Top (junction), Load at Bottom (ground connection)\n');
fprintf(fid, 'CE\n');

% ΓΕΩΜΕΤΡΙΑ
% GW 1: Οριζόντιος αγωγός
fprintf(fid, 'GW  1  %d  %.4f  0.0000  %.4f  %.4f  0.0000  %.4f  %.4f\n', ...
        segs_hor, -L_total/2, h, L_total/2, h, a_wire);

% GW 2: Κατακόρυφο σύρμα ΤΕΡΜΑΤΙΣΜΟΥ (Δεξιά)
fprintf(fid, 'GW  2  %d  %.4f  0.0000  %.4f  %.4f  0.0000  0.0000  %.4f\n', ...
        segs_vert, L_total/2, h, L_total/2, a_feed);

% GW 3: Κατακόρυφο σύρμα ΤΡΟΦΟΔΟΣΙΑΣ (Αριστερά)
fprintf(fid, 'GW  3  %d  %.4f  0.0000  %.4f  %.4f  0.0000  0.0000  %.4f\n', ...
        segs_vert, -L_total/2, h, -L_total/2, a_feed);

fprintf(fid, 'GE  1\n'); 
fprintf(fid, 'GN  1  0  0  0  0  0  0  0\n');   % Αρχικά Perfect Ground

%% ΑΝΤΙΣΤΑΣΗ ΤΕΡΜΑΤΙΣΜΟΥ (LOAD)
% Τοποθέτηση στο Tag 2, Segment 5 (Πάνω στο δεξί κατακόρυφο σύρμα, κοντά στη γη)
fprintf(fid, 'LD  0  2  5  5  %.1f  0  0\n', R_load);

%% ΠΗΓΗ ΔΙΕΓΕΡΣΗΣ ΚΑΙ ΣΥΧΝΟΤΗΤΑ
% Τοποθέτηση στο Tag 3, Segment 1 (Πάνω στο αριστερό κατακόρυφο σύρμα, στην αρχή του κυριώς αγωγού)
fprintf(fid, 'EX  0  3  1  0  1.0  0.0\n');

fprintf(fid, 'FR  0  1  0  0  %.2f  0\n', f0);
fprintf(fid, 'EN\n');
fclose(fid);

disp(['Το αρχείο ', filename, ' δημιουργήθηκε επιτυχώς!']);