% forward_stub.m
% theta = [R0, r0, pitch, cL_over_cT, j_amp]
% Return a struct with fields: m_pred, mu_pred, lamC_pred, Fres, Prad
function pred = forward_stub(theta)
    R0 = theta(1); 
    r0 = theta(2); 
    pitch = theta(3);
    cR = theta(4); 
    j = theta(5); %#ok<NASGU>

    alpha = 1/137.036;
    cT = 299792458; % or your c_T
    q  = 1.602176634e-19;

    % --- REPLACE everything below with your GDM computations ---
    lamC_pred = 2*pi*R0;
    I  = q * (cT/(2*pi*R0));
    mu_pred = I * (pi*R0^2);
    m_pred  = 9.0e-31 * ( (R0*1e12)^0.1 ) * ( (r0/R0)/alpha )^0.1 * (1+0.05*(cR-1.6068));
    Fres    = (r0/R0 - alpha) + 0.1*pitch;
    Prad    = 1e-15*(pitch^2);
    % -----------------------------------------------------------

    pred = struct('m_pred',m_pred,'mu_pred',mu_pred,'lamC_pred',lamC_pred,'Fres',Fres,'Prad',Prad);
end
