function initLm = getinitalLm(angle)
% start point muscle length
FCR_Lm_reg  = @(x) -0.0147*x + 0.0606;
FCU_Lm_reg  = @(x) -0.0124*x + 0.0412;
ECRL_Lm_reg = @(x) 0.0788 + 0.004631*cos(x*0.9933 ) +  0.0156*sin(x*0.9933 ) -0.0008986 *cos(2*x*0.9933 ) -0.003016*sin(2*x*0.9933);
ECRB_Lm_reg = @(x) 0.05642 + 0.00119*cos(x*0.9547) + 0.02014*sin(x*0.9547) -0.0002117 *cos(2*x*0.9547) -0.003375*sin(2*x*0.9547);
ECU_Lm_reg  = @(x) 0.05618 -0.00194*cos(x*1.286) + 0.005234*sin(x*1.286) +  0.0008253*cos(2*x*1.286) -0.0001549 *sin(2*x*1.286);

initLm = [FCR_Lm_reg(angle) FCU_Lm_reg(angle)  ECRL_Lm_reg(angle) ECRB_Lm_reg(angle) ECU_Lm_reg(angle)];
end














