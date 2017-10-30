%%% FBMC Implementation 
% To convert into System Generator

clear all;
close all;

%% System Parameters

numFFT = 256;
numGuards = 20;          % for both sides
K = 3;                   % Filter length
bitsPerSubCarrier = 4;   % 2: QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
numSymbols = 1000;       % Simulation length; OFDM Symbols

%% FBMC Filter

L = numFFT-2*numGuards; % No. of Carriers  % Number of complex symbols per OFDM symbol

HkOneSided = [0.911438 0.411438];
Hk = [fliplr(HkOneSided) 1 HkOneSided];

%BER def

KF = K*numFFT;
KL = K*L;

dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);
sumFBMCSpec = zeros(numFFT*K*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);

inpData = zeros( bitsPerSubCarrier*L/2, numSymbols); rxBits = inpData;
txSigAll = complex(zeros(KF, numSymbols));

%% FBMC modulation

for symbolNr = 1%:numSymbols
    
    % Half the needed symbols to account for oversampling by 2
    inpData(:, symbolNr) = randi([0 1], bitsPerSubCarrier*L/2, 1);
    modData = qam_mod_custom(inpData(:, symbolNr));
    
%     for Nidx = 1:2
        % OQAM Modulator: alternate real and imaginary
        if rem(symbolNr,2)==1 %Odd no of symbol
            dataSubCar(1:2:L) = real(modData);
            dataSubCar(2:2:L) = 1i*imag(modData);
        else %even no of symbol
            dataSubCar(1:2:L) = 1i*imag(modData);
            dataSubCar(2:2:L) = real(modData);
        end
        
        figure,subplot(2,1,1)
        plot(abs(dataSubCar));
        subplot(2,1,2)
        plot(angle(dataSubCar));
        title('datSubCar');
        
        % Upsample by K, pad with GB, and filter with the prototype filter
        dataSubCarUp(1:K:end) = dataSubCar;
        
        figure,subplot(2,1,1)
        plot(abs(dataSubCarUp));
        subplot(2,1,2)
        plot(angle(dataSubCarUp));
        title('dataSubCarUp');
        
        dataBitsUpPad = [zeros(numGuards*K,1); dataSubCarUp; zeros(numGuards*K,1)];
        X1 = filter(Hk, 1, dataBitsUpPad);

        % Remove 1/2 filter length delay
        X = [X1(K:end); zeros(K-1,1)];
        
        figure,subplot(2,1,1)
        plot(abs(X));
        subplot(2,1,2)
        plot(angle(X));
        title('X');

        % Compute IFFT of length KF => Transmitted signal
        scFactor = KF/sqrt(KL); 
        txSig = scFactor.*fftshift(ifft(X)); 
        % Combine the individual txSigs
%     end
    % Compute psd
%    step(specView, txSig);
    [specFBMC, fFBMC] = periodogram(txSig,rectwin(length(txSig)),numFFT*K*2,1);
    sumFBMCSpec = sumFBMCSpec + specFBMC;

    txSigAll(:,symbolNr) = txSig;
end

% Display Spectrums
sumFBMCSpec = sumFBMCSpec / mean(sumFBMCSpec(1+K+2*numGuards*K:end-2*numGuards*K-K));
figure,plot(fFBMC,10*log10(sumFBMCSpec)); hold on; grid on
axis([0 1 -130 10]);
xlabel('Normalized frequency'); ylabel('PSD (dBW/Hz)')
title(['FBMC with ' num2str(K) ' overlapping symbols'])
set(gcf, 'Position', [168 652 560 420]);


%% FBMC Receiver: no channel, SISO

for symbolNr = 1%:numSymbols
    rxSig = txSigAll(:, symbolNr);
    
    % Perform FFT
    Xrx = fft(fftshift(rxSig/scFactor));
    
    % Filter with prototype filter (matched filter)
    Xrx1 = filter(Hk, 1, Xrx);
    % Remove K-1 delay elements
    Xrx1 = [Xrx1(K:end); zeros(K-1,1)];
    % Remove guards
    X2 = Xrx1(numGuards*K+1:end-numGuards*K);
    
    % Downsample by 2K and extract real and imaginary part
    if rem(symbolNr, 2)
        % Imaginary part is K samples after real one
        r1 = real(X2(1:2*K:end));
        r2 = imag(X2(K+1:2*K:end));
        rcomb = complex(r1, r2)/K;
    else
        % Real part is K samples after imaginary one
        r1 = imag(X2(1:2*K:end));
        r2 = real(X2(K+1:2*K:end));
        rcomb = complex(r2, r1)/K;
    end
    
    % Perform hard decisions
    rxBits(:, symbolNr) = qam_demod_custom(rcomb);    
end
