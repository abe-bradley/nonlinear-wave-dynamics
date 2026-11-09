function v_out = noiser(v_in,SNR,scramble)
    %% Inputs
    % v_in and v_out are arrays containing complex numbers
    % SNR is the signal to noise ratio. It is defined as follows:
        % 1. Let v_in = sqrt(P0)*ones(1,N)
        % 2. Let P_out = abs(v_out).^2
        % 3. SNR = P0/std(P_out)

    %% Main
    Dv = size(v_in);

    % Apply amplitude noise
    v_out = v_in.*sqrt(1 + randn(Dv)/SNR);

    % Apply phase noise
    v_out = exp(2i*pi*randn(Dv)/SNR).*v_out;

    % Scramble phases
    if scramble
        v_out = exp(2i*pi*rand(Dv)).*v_out;
    end
end