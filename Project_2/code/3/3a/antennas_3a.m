% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775
% ΚΕΡΑΙΕΣ ΚΑΙ ΔΙΑΔΟΣΗ ΕΡΓΑΣΙΑ 2 - ΘΕΜΑ 3α

clc; clear; close all;

%% ΟΡΙΣΜΟΣ ΠΑΡΑΜΕΤΡΩΝ
f0 = 60;                    % Κεντρική συχνότητα (MHz)
lambda = 300 / f0;          % Μήκος κύματος 
L = 0.5 * lambda;           % Μήκος διπόλου 
a = (lambda / 200 ) / 2;    % Ακτίνα αγωγού 

% --- ΕΠΙΛΟΓΗ ΑΠΟΣΤΑΣΗΣ ---
% 1: S = lambda/100, 2: S = lambda/20, 3: S = lambda/4
case_idx = 1; % <--- Η ΑΛΛΑΓΗ ΓΙΝΕΤΑΙ ΕΔΩ (1, 2 ή 3)

if case_idx == 1
    S = lambda / 100;
    case_name = '1st_case_L100';
    segs_L = 31; 
    segs_S = 1;
elseif case_idx == 2
    S = lambda / 20;
    case_name = '2nd_case_L20';
    segs_L = 31; 
    segs_S = 3;
else
    S = lambda / 4;
    case_name = '3rd_case_L4';
    segs_L = 15; 
    segs_S = 7; % Αυξάνουμε εδώ για να πλησιάσει το μήκος των L segments
end


%% ΔΗΜΙΟΥΡΓΙΑ ΑΡΧΕΙΟΥ NEC
filename = sprintf('folded_dipole_geometry_%s.nec', case_name);
fid = fopen(filename, 'w');

fprintf(fid, 'CM Folded Dipole Antenna - f0=60MHz\n');
fprintf(fid, 'CM Designed by Karatis Dimitrios (10775)\n');
fprintf(fid, 'CE\n');

% --- Γεωμετρία (Άξονας Z) ---
% Αγωγός 1 (Τροφοδοσία) στο x = -S/2
fprintf(fid, 'GW  1  %d  %.4f  0  %.4f  %.4f  0  %.4f  %.4f\n', ...
        segs_L, -S/2, -L/2, -S/2, L/2, a);

% Αγωγός 2 (Παράλληλος) στο x = S/2
fprintf(fid, 'GW  2  %d  %.4f  0  %.4f  %.4f  0  %.4f  %.4f\n', ...
        segs_L, S/2, -L/2, S/2, L/2, a);

% Άνω βραχυκυκλωτήρας (z = L/2)
fprintf(fid, 'GW  3  %d  %.4f  0  %.4f  %.4f  0  %.4f  %.4f\n', ...
        segs_S, -S/2, L/2, S/2, L/2, a);

% Κάτω βραχυκυκλωτήρας (z = -L/2)
fprintf(fid, 'GW  4  %d  %.4f  0  %.4f  %.4f  0  %.4f  %.4f\n', ...
        segs_S, -S/2, -L/2, S/2, -L/2, a);

fprintf(fid, 'GE  0\n');

%% ΠΗΓΗ ΔΙΕΓΕΡΣΗΣ ΚΑΙ ΣΥΧΝΟΤΗΤΑ
% Πηγή στο κεντρικό segment του Tag 1
fprintf(fid, 'EX  0  1  %d  0  1.0  0.0\n', round(segs_L/2));

% Ορισμός της κεντρικής συχνότητας f0 = 60 MHz
fprintf(fid, 'FR  0  1  0  0  %.2f  0\n', f0);

fprintf(fid, 'EN\n');
fclose(fid);

disp(['Το αρχείο ', filename, ' δημιουργήθηκε επιτυχώς!']);