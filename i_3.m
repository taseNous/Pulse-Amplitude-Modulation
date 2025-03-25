clc;
clear;
close all;

%% Παράμετροι
M = 8; % Τάξη M-PAM (μπορεί επίσης να δοκιμαστεί με M=2)
log2M = log2(M); % Αριθμός bits ανά σύμβολο
Lb = 100000; % Συνολικός αριθμός bits
Lb = floor(Lb / log2M) * log2M; % Προσαρμογή ώστε να είναι διαιρετός με log2(M)
Es = 1; % Ενέργεια ανά σύμβολοl
Rs = 250e3; % Ρυθμός συμβόλων (250 Ksymbols/sec)
Ts = 1 / Rs; % Διάρκεια συμβόλου
fc = 2.5e6; % Συχνότητα φέρουσας (2.5 MHz)
Tc = 1 / fc; % Περίοδος φέρουσας
Tsamp = Tc / 4; % Περίοδος δειγματοληψίας (4 δείγματα ανά περίοδο φέρουσας)
Ns = round(Ts / Tsamp); % Αριθμός δειγμάτων ανά σύμβολο
SNR_dB_values = 0:2:20; % Εύρος SNR (dB)

%% Δημιουργία Δυαδικής Ακολουθίας Εισόδου
bits = randi([0 1], 1, Lb); % Τυχαία δυαδική ακολουθία

%% Mapping: Bits σε σύμβολα
symbols = reshape(bits, log2M, []); % Ομαδοποίηση bits σε log2(M) ανά σύμβολο
decSymbols = bi2de(symbols.', 'left-msb'); % Μετατροπή σε δεκαδικά σύμβολα
Am = (2 * (0:M-1) - (M - 1)) * sqrt(Es / (2 * (M - 1)^2)); % Επίπεδα πλάτους PAM

%% Διαμόρφωση παλμού
gT = ones(1, Ns); % Ορθογώνιος Παλμός
gT = gT / sqrt(sum(gT.^2)); % Κανονικοποίηση στη μονάδα ενέργειας

s_t = reshape(repmat(Am(decSymbols + 1), Ns, 1), 1, []); % Σήμα PAM με διαμόρφωση παλμού

%% Διαμόρφωση Σήματος
num_samples = length(s_t); % Συνολικός αριθμός δειγμάτων στο s_tt
t_carrier = linspace(0, (num_samples-1)*Tsamp, num_samples); % Χρονικό διάνυσμα για τη φέρουσα
carrier = cos(2 * pi * fc * t_carrier); % Κυματομορφή φέρουσας
transmitted_signal = s_t .* carrier; % Διαμορφωμένο σήμα

%% Οπτικοποίηση Μεταδιδόμενων Συμβόλων (Αντιστοίχιση)
figure;
stem(1:20, decSymbols(1:20), 'filled', 'LineWidth', 1.5);
title('Transmitted Symbols (Mapping)');
xlabel('Symbol Index');
ylabel('Amplitude');
grid on;

%% Αρχικοποίηση BER και SER
BER_values = zeros(size(SNR_dB_values)); % Αποθήκευση BER
SER_values = zeros(size(SNR_dB_values)); % Αποθήκευση SER

%% Υπολογισμός BER και SER
for idx = 1:length(SNR_dB_values)
    SNR_dB = SNR_dB_values(idx);
    SNR_linear = 10^(SNR_dB / 10);
    
    % Προσθήκη AWGN
    No = Es / (SNR_linear * log2M); % Ισχύς φασματικής πυκνότητας θορύβου
    sigma = sqrt(No / 2); % Τυπική απόκλιση θορύβου
    noise = sigma * randn(size(transmitted_signal));
    received_signal = transmitted_signal + noise;
    
    % Αποδιαμόρφωση
    demodulated_signal = received_signal .* carrier; % Πολλαπλασιασμός με φέρουσα
    
    % Ολοκλήρωση
    integrated_signal = zeros(1, length(decSymbols)); % Αρχικοποίηση ολοκληρωμένου σήματος
    for i = 1:length(decSymbols)
        start_idx = (i-1)*Ns + 1;
        end_idx = i*Ns;
        integrated_signal(i) = sum(demodulated_signal(start_idx:end_idx) .* gT);
    end
    
    % Κλιμάκωση ολοκληρωμένου σήματος
    integrated_signal = integrated_signal / max(abs(integrated_signal)) * max(abs(Am));
    
    % Ανίχνευση Συμβόλων
    detected_symbols = zeros(size(integrated_signal));
    for i = 1:length(integrated_signal)
        [~, detected_symbols(i)] = min(abs(integrated_signal(i) - Am)); % Ανίχνευση πλησιέστερου πλάτους
    end
    detected_symbols = detected_symbols - 1; % Προσαρμογή στο εύρος [0, M-1]
    
    % Οπτικοποίηση Ληφθέντων Συμβόλων (Αντιστοίχιση) για SNR = 20 dB
    if SNR_dB == 20
        figure;
        stem(1:20, detected_symbols(1:20), 'filled', 'LineWidth', 1.5);
        title('Ληφθέντα Σύμβολα (Αντιστοίχιση)');
        xlabel('Δείκτης Συμβόλου');
        ylabel('Πλάτος');
        grid on;
    end
    
    % Αντιστοίχιση: Σύμβολα σε Bits
    detected_bits = de2bi(detected_symbols, log2M, 'left-msb')'; % Μετατροπή συμβόλων σε bits
    detected_bits = detected_bits(:)'; % Επίπεδη αναπαράσταση σε 1D πίνακα
    
    % Υπολογισμός BER
    if ~isempty(detected_bits)
        minLength = min(length(bits), length(detected_bits)); % Ευθυγράμμιση μηκών
        BER_values(idx) = sum(bits(1:minLength) ~= detected_bits(1:minLength)) / Lb; % Κανονικοποίηση με τον συνολικό αριθμό bits
    else
        BER_values(idx) = NaN; % Ανάθεση NaN αν η ανίχνευση αποτύχει
        fprintf('Η ανίχνευση απέτυχε για SNR = %d dB\n', SNR_dB);
    end
    
    % Υπολογισμός SER
    decSymbols = decSymbols(:)'; % Διασφάλιση ότι είναι διάνυσμα γραμμής
    detected_symbols = detected_symbols(:)'; % Διασφάλιση ότι είναι διάνυσμα γραμμής
    num_symbol_errors = sum(detected_symbols ~= decSymbols); % Μέτρηση λαθών συμβόλων
    SER_values(idx) = num_symbol_errors / length(decSymbols); % Αναλογία λανθασμένων συμβόλων
end

%% Debugging
disp('First 10 transmitted symbols:');
disp(decSymbols(1:10));

disp('First 10 detected symbols (SNR = 20 dB):');
disp(detected_symbols(1:10));

fprintf('Amplitude levels: [%f, %f]\n', min(Am), max(Am));
fprintf('Integrated signal range: [%f, %f]\n', min(integrated_signal), max(integrated_signal));

figure;
semilogy(SNR_dB_values, BER_values, '-o', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB_values, SER_values, '-s', 'LineWidth', 1.5);
title('BER and SER vs. SNR');
xlabel('SNR (dB)');
ylabel('Error Rate');
legend('BER', 'SER');
grid on;