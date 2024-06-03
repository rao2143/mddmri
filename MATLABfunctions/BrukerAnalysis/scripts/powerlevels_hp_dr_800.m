clear all
%clc

%Accept specs 3.2 mm
p901H = 2.44; plw1H = 216.4;
p9013C = 3.38; plw13C = 170.0;
p9031P = 3.03; plw31P = 200;
p9079Br = 3.55; plw79Br = 151.1;

B1ref1H = 204/2;             % kHz
plwref1H = 200;                     % Watt
B1ref13C = 186/2;           % kHz
plref13C = 200;                   % Watt
B1ref31P = 191/2;           % kHz
plref31P = 200;                   % Watt

spinrate = 18e3;

B1hp1H = 100;   % kHz
B1dec = 62;                   % kHz
B1hp13C = B1hp1H;   % kHz
B1hp31P = B1hp1H;   % kHz
B1cp13C = 60;          % kHz
B1cp31P = 60;          % kHz
B1spinlock13C = 52;          % kHz

p90hp = 1000/B1hp1H/4                    % microsecs
p180dec = 1000/B1dec/2; %us
pcpd3 = p180dec
p180rec = 1/(spinrate*9*2)*10^6 % microsecs
B1rec = 1/(2*p180rec*0.001);   % kHz

plw1 = plref13C*(B1hp13C/B1ref13C)^2   % Watt
plw21 = plref13C*(B1cp13C/B1ref13C)^2  % Watt
plw2 = plref31P*(B1hp31P/B1ref31P)^2   % Watt
plw22 = plref31P*(B1cp31P/B1ref31P)^2  % Watt

plw3 = plwref1H*(B1hp1H/B1ref1H)^2                % Watt
plw12 = plwref1H*(B1dec/B1ref1H)^2                % Watt
plw32 = plwref1H*(B1rec/B1ref1H)^2                % Watt
spw0 = plw3
spw2 = plw3

plw13 = plref13C*(B1spinlock13C/B1ref13C)^2  % Watt
in_f = 9*4*p180rec

