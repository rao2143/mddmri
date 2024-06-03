clear all
%clc

B1ref1H = 148/2;             % kHz
plwref1H = 100;                     % Watt
B1ref13C = 178/2;           % kHz
plref13C = 250;                   % Watt
B1ref31P = 216/2;           % kHz
plref31P = 250;                   % Watt

spinrate = 12e3;

B1hp1H = 100;   % kHz
B1dec = 62;                   % kHz
B1hp13C = B1hp1H;   % kHz
B1hp31P = B1hp1H;   % kHz
B1cp13C = 70;          % kHz
B1cp31P = 70;          % kHz

p90hp = 1000/B1hp1H/4                    % microsecs
p180dec = 1000/B1dec/2; %us
pcpd3 = p180dec
%pl180dec = 7.35;                 % microsec
%B1dec = 1/(2*pl180dec*0.001)   % kHz

plw1 = plref13C*(B1hp13C/B1ref13C)^2   % Watt
plw21 = plref13C*(B1cp13C/B1ref13C)^2  % Watt
plw2 = plref31P*(B1hp31P/B1ref31P)^2   % Watt
plw22 = plref31P*(B1cp31P/B1ref31P)^2  % Watt

plw3 = plwref1H*(B1hp1H/B1ref1H)^2                % Watt
plw12 = plwref1H*(B1dec/B1ref1H)^2                % Watt
spw0 = plw3
spw2 = plw3

