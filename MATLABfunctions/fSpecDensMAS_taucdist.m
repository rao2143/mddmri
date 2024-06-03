function J = fSpecDensMAS_taucdist(omega,omegaR,tauc,w)

if isstruct(omega)
    fields = fieldnames(omega);
    
    for nfield = 1:numel(fields)
        field = fields{nfield};
        J.(field) =fSpecDensMAS(omega.(field),omegaR,tauc')*w;
    end
    
else
    J =fSpecDensMAS(omega,omegaR,tauc')*w;
end
        
