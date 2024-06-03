function Y = fSumLorentz(Pin,Xin,Pnorm,Xnorm,Ynorm)if nargin == 2    Pnorm = 1;    Xnorm = 1;    Ynorm = 1;endPin = Pin.*Pnorm;nu = Xin*Xnorm;npeak = length(Pin)/3;T2 = Pin(1:npeak);ampl = Pin(npeak+1:2*npeak);offset = Pin(2*npeak+1:3*npeak);phase = 0;Icalc = zeros(size(nu));for npeak = 1:length(ampl)	Icalc = Icalc + fLorentz(nu,ampl(npeak),T2(npeak),offset(npeak),phase);endIcalc = real(Icalc);Y = Icalc/Ynorm;