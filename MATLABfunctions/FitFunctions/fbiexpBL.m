function Y = fbiexpBL(Pin,Xin,Pnorm,Xnorm,Ynorm);Pin = Pin.*Pnorm;Xin = Xin*Xnorm;Y0 = Pin(1);Df = Pin(2);Ds = Pin(3);ff = Pin(4);BL = Pin(5);Y = Y0.*(ff.*exp( - Df.*Xin) + (1-ff).*exp( - Ds.*Xin)) + BL;Y = Y/Ynorm;