function Y = fLNs(Pin,Xin,Pnorm,Xnorm,Ynorm);Pin = Pin.*Pnorm;Xin = Xin*Xnorm;I0 = Pin(1);Dmean = Pin(2);sigma = Pin(3);Dspike = Pin(4);fspike = Pin(5);D0 = Dmean/exp(sigma^2/2);Dmax = D0/exp(sigma^2);N = 100;D = logspace(log10(Dmax)-3*sigma,log10(Dmax)+3*sigma,N);[karray, Darray] = ndgrid(Xin,D);SD = 1./(Darray.*sigma*sqrt(2*pi)).*exp(-1/2*((log(Darray)-log(D0))./sigma).^2).*exp(-karray.*Darray);mSD = (SD(:,1:(N-1)) + SD(:,2:N))/2;dD = D(2:N) - D(1:(N-1));Y = I0.*((1-fspike).*mSD*dD' + fspike.*exp(-Xin.*Dspike));Y = Y/Ynorm;