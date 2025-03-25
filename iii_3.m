clc;
clear;
close all;

%% Παράμετροι
M_values = [2, 8]; % Επίπεδα M-PAM
log2M_values = log2(M_values); % Αριθμός bits ανά σύμβολο για κάθε M
Lb = 100000; % Συνολικός αριθμός bits
Es = 1; % Ενέργεια ανά σύμβολο
Rs = 250e3; % Ρυθμός συμβόλων (250 Ksymbols/sec)
Ts = 1 / Rs; % Διάρκεια συμβόλου
fc = 2.5e6; % Συχνότητα φέρουσας (2.5 MHz)
Tc = 1 / fc; % Περίοδος φέρουσας
Tsamp = Tc / 4; % Περίοδος δειγματοληψίας (4 δείγματα ανά περίοδο φέρουσας)
SNR_dB_values = 0:2:20; % Εύρος SNR (dB)

%% Αποθήκευση SER
SER_values = zeros(length(SNR_dB_values), length(M_values)); % Αποθήκευση SER για κάθε M

%% Επανάληψη για τις τιμές M
for m_idx = 1:length(M_values)
    M = M_values(m_idx);
    log2M = log2M_values(m_idx);
    Lb = floor(Lb / log2M) * log2M; % Προσαρμογή ώστε να είναι διαιρετός με log2(M)
    bits = randi([0 1], 1, Lb); % Τυχαία δυαδική ακολουθία

    % Αντιστοίχιση: Bits σε Σύμβολα
    symbols = reshape(bits, log2M, []); % Ομαδοποίηση bits σε log2(M) ανά σύμβολο
    decSymbols = bi2de(symbols.', 'left-msb'); % Μετατροπή σε δεκαδικά σύμβολα
    Am = (2 * (0:M-1) - (M - 1)) * sqrt(Es / (2 * (M - 1)^2)); % Επίπεδα πλάτους PAM
    
    % Διαμόρφωση Παλμού
    Ns = round(Ts / Tsamp); % Αριθμός δειγμάτων ανά σύμβολο
    gT = ones(1, Ns); % Ορθογώνιος παλμός
    gT = gT / sqrt(sum(gT.^2)); % Κανονικοποίηση σε μοναδιαία ενέργεια
    s_t = reshape(repmat(Am(decSymbols + 1), Ns, 1), 1, []); % Σήμα PAM με διαμόρφωση παλμού
    
    % Διαμόρφωση
    t_carrier = linspace(0, (length(s_t)-1)*Tsamp, length(s_t)); % Χρονικό διάνυσμα για τη φέρουσα
    carrier = cos(2 * pi * fc * t_carrier); % Κυματομορφή φέρουσας
    transmitted_signal = s_t .* carrier; % Διαμορφωμένο σήμα
    
    %% Υπολογισμός SER για κάθε τιμή SNR
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
        integrated_signal = integrated_signal / max(abs(integrated_signal)) * max(abs(Am)); % Κλιμάκωση σήματος
    
        % Ανίχνευση Συμβόλων
        detected_symbols = zeros(size(integrated_signal));
        for i = 1:length(integrated_signal)
            [~, detected_symbols(i)] = min(abs(integrated_signal(i) - Am)); % Ανίχνευση πλησιέστερου πλάτους
        end
        detected_symbols = detected_symbols - 1; % Προσαρμογή στο εύρος [0, M-1]
    
        % Διασφάλιση ευθυγράμμισης διανυσμάτων
        decSymbols = decSymbols(:)'; % Μετατροπή σε διάνυσμα γραμμής
        detected_symbols = detected_symbols(:)'; % Μετατροπή σε διάνυσμα γραμμής
        assert(length(decSymbols) == length(detected_symbols), 'Ασυμφωνία στα μήκη διανυσμάτων συμβόλων.');
    
        % Υπολογισμός αριθμού λαθών συμβόλων
        num_symbol_errors = sum(detected_symbols ~= decSymbols); % Μέτρηση λαθών συμβόλων
    
        % Υπολογισμός SER
        SER_values(idx, m_idx) = num_symbol_errors / length(decSymbols); % Ποσοστό λανθασμένων συμβόλων
    end
end

%% Plot SER
figure;
hold on;
for m_idx = 1:length(M_values)
    semilogy(SNR_dB_values, SER_values(:, m_idx), '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('%d-PAM', M_values(m_idx)));
end
title('SER vs. SNR for M = 2 and M = 8 (Simple Mapping)');
xlabel('SNR (dB)');
ylabel('Symbol Error Rate (SER)');
legend show;
grid on;
ylim([1e-5, 1]);
