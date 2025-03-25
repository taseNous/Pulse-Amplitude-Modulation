clc;
clear;
close all;

%% Παράμετροι
M_values = [2, 8]; % Επίπεδα M-PAM
log2M_values = log2(M_values); % Αριθμός bits ανά σύμβολο για κάθε M
Lb = 100000; % Συνολικός αριθμός bits
Es = 1; % Αριθμός bits ανά σύμβολο για κάθε M
Rs = 250e3; % Ρυθμός συμβόλων (250 Ksymbols/sec)
Ts = 1 / Rs; % Διάρκεια συμβόλου
fc = 2.5e6; % Συχνότητα φέρουσας (2.5 MHz)
Tc = 1 / fc; % Περίοδος φέρουσας
Tsamp = Tc / 4; % Περίοδος δειγματοληψίας (4 δείγματα ανά περίοδο φέρουσας
Ns = round(Ts / Tsamp); % Δείγματα ανά σύμβολο
SNR_dB_values = 0:2:20; % Εύρος SNR (dB)

%% Αρχικοποίηση αποθήκευσης BER
BER_simple = zeros(length(SNR_dB_values), length(M_values));
BER_gray = zeros(length(SNR_dB_values), 1); % Μόνο για M = 8 με κωδικοποίηση Gray

%% Προσομοίωση για κάθε M
for m_idx = 1:length(M_values)
    M = M_values(m_idx);
    log2M = log2M_values(m_idx);
    Lb_adj = floor(Lb / log2M) * log2M; % Προσαρμογή αριθμού bits
    bits = randi([0 1], 1, Lb_adj); % Τυχαία δυαδική ακολουθία

    % Αντιστοίχιση Συμβόλων (Απλή)
    symbols = reshape(bits, log2M, []); % Ομαδοποίηση bits για κάθε σύμβολο
    decSymbols_simple = bi2de(symbols.', 'left-msb'); % Δεκαδικά σύμβολα (Απλή)
    Am_simple = (2 * (0:M-1) - (M-1)) * sqrt(Es / (2 * (M-1)^2));
    gray_map = bitxor(0:M-1, floor((0:M-1)/2)); % Δημιουργία αντιστοίχισης κώδικα Gray
    Am_gray = Am_simple(gray_map + 1); % Αντιστοίχιση πλατών χρησιμοποιώντας τον κώδικα Gray
    
    % Αντιστοίχιση Συμβόλων (Gray για M = 8)
    if M == 8
        gray_map = bitxor(0:M-1, floor((0:M-1)/2)); % Αντιστοίχιση κώδικα Gray
        Am_gray = Am_simple; % Ίδια επίπεδα πλάτους
        graySymbols = gray_map(decSymbols_simple + 1); % Εφαρμογή αντιστοίχισης Gray
    end
    
    % Διαμόρφωση PAM (Διαμόρφωση παλμού)
    gT = ones(1, Ns) / sqrt(Ns); % Κανονικοποιημένος ορθογώνιος παλμός
    s_t_simple = reshape(repmat(Am_simple(decSymbols_simple + 1), Ns, 1), 1, []); % Απλή αντιστοίχιση
    if M == 8
        s_t_gray = reshape(repmat(Am_gray(graySymbols + 1), Ns, 1), 1, []); % Αντιστοίχιση Gray
    end
    
    % Διαμόρφωση φέρουσας
    num_samples = length(s_t_simple);
    t_carrier = linspace(0, (num_samples-1)*Tsamp, num_samples); % Χρονικό διάνυσμα
    carrier = cos(2 * pi * fc * t_carrier); % Κυματομορφή φέρουσας
    transmitted_signal_simple = s_t_simple .* carrier;
    if M == 8
        transmitted_signal_gray = s_t_gray .* carrier;
    end
    
    % Υπολογισμός BER
    for snr_idx = 1:length(SNR_dB_values)
        SNR_dB = SNR_dB_values(snr_idx);
        SNR_linear = 10^(SNR_dB / 10);
    
        % Προσθήκη θορύβου
        No = Es / (SNR_linear * log2M); % Ισχύς φασματικής πυκνότητας θορύβου
        sigma = sqrt(No / 2); % Τυπική απόκλιση θορύβου AWGN
        noise_simple = sigma * randn(size(transmitted_signal_simple)); % Θόρυβος για απλή αντιστοίχιση
        if M == 8
            noise_gray = sigma * randn(size(transmitted_signal_gray)); % Θόρυβος για αντιστοίχιση Gray
        end
    
        received_signal_simple = transmitted_signal_simple + noise_simple;
        if M == 8
            received_signal_gray = transmitted_signal_gray + noise_gray;
        end
    
        % Αποδιαμόρφωση και Ανίχνευση (Απλή)
        demodulated_simple = received_signal_simple .* carrier; % Πολλαπλασιασμός με φέρουσα
        integrated_simple = zeros(1, length(decSymbols_simple));
        for i = 1:length(decSymbols_simple)
            start_idx = (i-1)*Ns + 1;
            end_idx = i*Ns;
            integrated_simple(i) = sum(demodulated_simple(start_idx:end_idx) .* gT);
        end
        % Ανίχνευση πλησιέστερου επιπέδου πλάτους
        integrated_simple = integrated_simple / max(abs(integrated_simple)) * max(abs(Am_simple));
        [~, detected_simple] = min(abs(integrated_simple(:) - Am_simple), [], 2); % Διασφάλιση ότι το Am_simple είναι διάνυσμα γραμμής
        detected_simple = detected_simple - 1; % Προσαρμογή στο εύρος [0, M-1]
        
        % Μετατροπή ανιχνευμένων συμβόλων σε δυαδικά bits
        detected_bits = de2bi(detected_simple, log2M, 'left-msb')'; % Μετατροπή σε δυαδικό
        detected_bits = detected_bits(:)'; % Μετατροπή σε διάνυσμα γραμμής
        
        % Ευθυγράμμιση μηκών και υπολογισμός BER
        minLength = min(length(bits), length(detected_bits)); % Ευθυγράμμιση μηκών
        BER_simple(snr_idx, m_idx) = sum(bits(1:minLength) ~= detected_bits(1:minLength)) / minLength;
    
        detected_simple = detected_simple - 1; % Προσαρμογή στο εύρος [0, M-1]
    
        % BER για Απλή Αντιστοίχιση
        detected_bits = detected_bits(:)'; % Μετατροπή ανιχνευμένων bits σε διάνυσμα γραμμής
        minLength = min(length(bits), length(detected_bits)); % Ευθυγράμμιση μηκών
        BER = sum(bits(1:minLength) ~= detected_bits(1:minLength)) / Lb;
    
    
        % Αποδιαμόρφωση και Ανίχνευση (Gray, μόνο για M = 8)
        if M == 8
            demodulated_gray = received_signal_gray .* carrier;
            integrated_gray = zeros(1, length(graySymbols));
            for i = 1:length(graySymbols)
                start_idx = (i-1)*Ns + 1;
                end_idx = i*Ns;
                integrated_gray(i) = sum(demodulated_gray(start_idx:end_idx) .* gT);
            end
            integrated_gray = integrated_gray / max(abs(integrated_gray)) * max(abs(Am_gray));
            [~, detected_gray] = min(abs(integrated_gray(:) - Am_gray), [], 2);
            detected_gray = detected_gray - 1; % Προσαρμογή στο [0, M-1]
            % Αντιστροφή αντιστοίχισης Gray
            inv_gray_map = zeros(1, M);
            inv_gray_map(gray_map + 1) = 0:M-1;
            decoded_gray = inv_gray_map(detected_gray + 1); % Χαρτογράφηση ανιχνευμένων συμβόλων πίσω στα αρχικά
    
            detected_bits_gray = de2bi(decoded_gray, log2M, 'left-msb')';
            detected_bits_gray = detected_bits_gray(:)';
            minLength_gray = min(length(bits), length(detected_bits_gray)); % Ευθυγράμμιση μηκών
            BER_gray(snr_idx) = sum(bits(1:minLength_gray) ~= detected_bits_gray(1:minLength_gray)) / Lb_adj;
    
        end
    end
    
    disp('BER για M=2:');
    disp(BER_simple(:, 1));

end

%% Plot
figure;
semilogy(SNR_dB_values, BER_simple(:, 1), '-o', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB_values, BER_simple(:, 2), '-s', 'LineWidth', 1.5);
semilogy(SNR_dB_values, BER_gray, '-d', 'LineWidth', 1.5);
title('BER vs. SNR for M = 2 & 8');
xlabel('SNR (dB)');
ylabel('BER');
legend('M=2 (Simple)', 'M=8 (Simple)', 'M=8 (Gray)');
grid on;