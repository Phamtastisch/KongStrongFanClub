clc;
clear all;
close all;



MUT_Luft = MUT(sparameters('Luftmessung.s2p'));
MUT_destWasser = MUT(sparameters('destWasser.s2p'));
MUT_kochsalzLoesung = MUT(sparameters('kochsalzLösung.s2p'));



% =========================================================
%  3-Punkt-Kalibrierung – Rational Function Model (RFM)
%  Referenzen: Luft, Wasser, Kochsalzlösung
% =========================================================

%% --- 1. Frequenzvektor (in Hz) ---

f= MUT_Luft.sobj_struct.Frequencies;

%% --- 2. Bekannte Referenz-Permittivitäten ---

% Luft: real, konstant
eps_Luft = ones(size(f));

% Wasser (Single-Debye, 25 °C)

% ToDo: das für Destilliertes Wasser herausfinden.

eps_s_w   = 78.4;      % statische Permittivität
eps_inf_w =  5.2;      % Hochfrequenz-Grenzwert
tau_w     = 9.4e-12;   % Relaxationszeit [s]
eps_water = eps_inf_w + (eps_s_w - eps_inf_w) ./ (1 + 1j*2*pi*f*tau_w);

% Kochsalzlösung (0.9 % NaCl, 25 °C)
% → Debye-Term + ionischer Leitfähigkeitsterm
eps_s_s   = 76.5;
eps_inf_s =  4.9;
tau_s     = 9.0e-12;
sigma_s   =  1.5;      % ionische Leitfähigkeit [S/m] — abhängig von Konzentration!
eps0      = 8.854e-12;
eps_kochsalzLoesung = eps_inf_s + (eps_s_s - eps_inf_s) ./ (1 + 1j*2*pi*f*tau_s) ...
    - 1j * sigma_s ./ (2*pi*f*eps0);

%% --- 3. Gemessene S11 der Referenzen (komplex) ---
% S11_air, S11_water, S11_saline: Spaltenvektoren über f
% → aus deinem Import bereits vorhanden

S11_Luft = MUT_Luft.S11_param;
S11_destWasser = MUT_destWasser.S11_param;
S11_kochsalzLoesung = MUT_kochsalzLoesung.S11_param;


%% --- 4. Kalibrierkoeffizienten A, B, D bei jeder Frequenz lösen ---
N = length(f);
A_cal = zeros(N,1);
B_cal = zeros(N,1);
D_cal = zeros(N,1);

for k = 1:N
    % Komplexe Logarithmen der Referenz-S11
    L1 = log(S11_Luft(k));
    L2 = log(S11_destWasser(k));
    L3 = log(S11_kochsalzLoesung(k));

    e1 = eps_Luft(k);
    e2 = eps_water(k);
    e3 = eps_kochsalzLoesung(k);

    % Lineares System:  M * [A; B; D] = rhs
    %   Herleitung: eps_i*(L_i + D) = A*L_i + B
    %   → A*L_i - eps_i*D + B = eps_i*L_i
    M = [ L1,  1,  -e1 ;
        L2,  1,  -e2 ;
        L3,  1,  -e3 ];

    rhs = [ e1*L1 ; e2*L2 ; e3*L3 ];

    % Konditionszahl prüfen (optional, aber empfehlenswert)
    if cond(M) > 1e10
        warning('Schlecht konditioniertes System bei f = %.3f GHz', f(k)/1e9);
    end

    x = M \ rhs;
    A_cal(k) = x(1);
    B_cal(k) = x(2);
    D_cal(k) = x(3);
end

%% --- 5. Permittivität der Probe berechnen ---
% S11_DUT: komplexer S11-Vektor der unbekannten Probe

MUT_Kongstrong_mitZucker = MUT(sparameters('KongStrong_mitZucker.s2p'))
MUT_Kongstrong_ohneZucker = MUT(sparameters('KongStrong_ohneZucker.s2p'))
MUT_Kongstrong_mitGeschmack = MUT(sparameters('KongStrong_mitGeschmack.s2p'))

S11_Kongstrong_mitZucker = MUT_Kongstrong_mitZucker.S11_param;
S11_Kongstrong_ohneZucker = MUT_Kongstrong_ohneZucker.S11_param;
S11_Kongstrong_mitGeschmack = MUT_Kongstrong_mitGeschmack.S11_param;
%Todo: Funktion schreiben, die MUT

eps_Kongstrong_mitZucker = eps_from_MUT(A_cal,B_cal,D_cal, ... 
                           S11_Kongstrong_mitZucker,N)
eps_Kongstrong_ohneZucker = eps_from_MUT(A_cal,B_cal,D_cal, ... 
    S11_Kongstrong_ohneZucker,N)
eps_Kongstrong_mitGeschmack = eps_from_MUT(A_cal,B_cal,D_cal, ... 
    S11_Kongstrong_mitGeschmack,N)

eps_array = {

    eps_Kongstrong_ohneZucker, 'eps Kongstrong ohne Zucker' , 'g'
    eps_Kongstrong_mitZucker,  'eps Kongstrong mit Zucker' , 'r'
    eps_Kongstrong_mitGeschmack, 'eps Kongstrong mit Geschmack' , 'y'

   
}

% %% --- 6. Plotten ---
% figure; hold on;
% 
% for i = 1:size(eps_array,1) 
% 
% subplot(2,1,1)
% plot(f/1e9, real(eps_array{i,1}), ...
%     'Color', eps_array{i,3}, ...
%     'DisplayName', eps_array{i,2},...
%     LineWidth=1.5); grid on
% 
% 
% xlabel('Frequenz (GHz)'); ylabel("\epsilon_r'")
% title('Realteil der Permittivität (Probe)')
% 
% subplot(2,1,2)
% plot(f/1e9, -imag(eps_array{i,1}), ...
%     'Color', eps_array{i,3}, ...
%     'DisplayName', eps_array{i,2},...
%     LineWidth=1.5); grid on
% 
% xlabel('Frequenz (GHz)'); ylabel("\epsilon_r''")
% title('Imaginärteil / Verlustfaktor')
% legend show
% end

%% --- 6. Plotten ---
figure;

% Erzeuge beide Subplots einmal und aktiviere hold für jede Achse
ax1 = subplot(2,1,1);
hold(ax1, 'on');
xlabel('Frequenz (GHz)'); ylabel("\epsilon_r'")
title('Realteil der Permittivität (Probe)')
grid(ax1, 'on');

ax2 = subplot(2,1,2);
hold(ax2, 'on');
xlabel('Frequenz (GHz)'); ylabel("\epsilon_r''")
title('Imaginärteil / Verlustfaktor')
grid(ax2, 'on');

% Plot in der Schleife (verwende geschweifte Klammern für Zellinhalt)
for i = 1:size(eps_array,1)
    plot(ax1, f/1e9, real(eps_array{i,1}), ...
        'Color', eps_array{i,3}, ...
        'DisplayName', eps_array{i,2}, ...
        'LineStyle','--',...
        'LineWidth', 1);

    plot(ax2, f/1e9, -imag(eps_array{i,1}), ...
        'Color', eps_array{i,3}, ...
        'DisplayName', eps_array{i,2}, ...
        'LineStyle','--',...
        'LineWidth', 1);
end

% Legende einmalig hinzufügen
legend(ax1, 'show', 'Location', 'best');
legend(ax2, 'show', 'Location', 'best');

