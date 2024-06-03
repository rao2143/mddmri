clear all

%Efree
B1ref1H = 148/2;             % kHz
plwref1H = 100;                     % Watt
B1ref13C = 181/2;           % kHz
plref13C = 250;                   % Watt

%CP-MAS
B1ref1H = 174/2;             % kHz
plwref1H = 100;                     % Watt
B1ref13C = 138/2;           % kHz
plref13C = 250;                   % Watt

spinrate = 5e3;

B1hp1H = 80e-3;   % kHz
B1dec = 62;                   % kHz
B1rec = 45;
B1hp13C = B1hp1H;   % kHz
B1cp13C = 60;          % kHz
B1spinlock13C = 24;   B1spinlock13C = linspace(24,48,6) % kHz

p90hp = 1000/B1hp1H/4                    % microsecs
p180dec = 1000/B1dec/2; %us
pcpd2 = p180dec
p180rec = 1/(spinrate*9*2)*1e6 % microsecs

plw1 = plref13C*(B1hp13C/B1ref13C)^2   % Watt
plw21 = plref13C*(B1cp13C/B1ref13C)^2  % Watt
plw13 = plref13C*(B1spinlock13C/B1ref13C).^2  % Watt

plw2 = plwref1H*(B1hp1H/B1ref1H)^2                % Watt
plw12 = plwref1H*(B1dec/B1ref1H)^2                % Watt
plw32 = plwref1H*(B1rec/B1ref1H)^2                % Watt
spw0 = plw2 % 1H cp same as hp

