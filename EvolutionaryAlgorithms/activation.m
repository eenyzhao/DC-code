% activaion dynamics
function activation = activation(afactor,e,threshold)

% emd_shift = zeros(size(e,1),1);
% B = zeros(sampleLength,5);
% emd = 6; % shift units -- represent the elecmechanical delay

% for zz = 1:5
%     B(emd+1:end,zz) = e(1:end-emd,zz);
% end   
for col = 1:size(e,2)
    % Activation signal 
    a(:,col) = (( exp (afactor* e(:,col)) - 1)./(exp(afactor)- 1 ) );
end

% thresholding
for col = 1:size(a,2)
    for row = 1:size(a,1)
        if  a(row,col) < threshold
            a(row,col) = 0; 
    end
end

activation = zeros(length(a),5);
%% moving average with 50ms 
% for k = 1:size(a,2)
%     activation(:,k) = movmean(a(:,k),50);
% end    
    
activation = a;
end 