newDP = dp';
nrows = size(newDP, 1);
widths = zeros(nrows,4);
nonlinear_energies = zeros(nrows,4);
peak_start = zeros(nrows,4);
peak_end = zeros(nrows,4);
delay_width = zeros(nrows,4);
delay_start = zeros(nrows,4);
delay_end = zeros(nrows,4);
valley_width = zeros(nrows,4);
valley_start = zeros(nrows,4);
valley_end = zeros(nrows,4);

warning off dg_spikewidth:novalley
warning off dg_spikewidth:novalley2
warning off dg_spikewidth:earlypeak

for wire = 1:4
    disp('NEW WIRE');
    progress = '';
    startIdx = (wire - 1)*32 + 1;
    endIdx = 32*wire;
    wire_data = newDP(:, startIdx:endIdx);
    for i=1:nrows
        if mod(i,10000)==0
            if strcmp(progress,'..........')==1
                progress = '.';
            else
                progress = strcat(progress,'.');
            end
            disp(progress);
        end
        [peak_width_vals, delay_vals, valley_vals] = dg_spikewidth(wire_data(i,:), 'allpoints');
        widths(i,wire) = peak_width_vals(1);
        t_0 = round(peak_width_vals(2));
        t_1 = round(peak_width_vals(3));
        
        peak_start(i,wire) = t_0;
        peak_end(i,wire) = t_1;
        delay_width(i,wire) = delay_vals(1);
        delay_start(i,wire) = delay_vals(2);
        delay_end(i,wire) = delay_vals(3);
        valley_width(i,wire) = valley_vals(1);
        valley_start(i,wire) = valley_vals(2);
        valley_end(i,wire) = valley_vals(3);
        
        if isnan(widths(i,wire)) == 0
            squared_voltage = wire_data(i,t_0:t_1).*wire_data(i,t_0:t_1);
            if t_0 ~= 1 && t_1 ~= 32
                nonlinear_energies(i,wire) = sum(squared_voltage...
                    - wire_data(i,t_0-1:t_1-1).*wire_data(i,t_0+1:t_1+1))/widths(i,wire);
            elseif t_1 ~= 32
                nonlinear_energies(i,wire) = sum(squared_voltage...
                    - [wire_data(i,t_0), wire_data(i,t_0:t_1-1)].*wire_data(i,t_0+1:t_1 + 1))/widths(i,wire);
            elseif t_0 ~= 1
                nonlinear_energies(i,wire) = sum(squared_voltage...
                    - wire_data(i,t_0-1:t_1-1).*[wire_data(i,t_0+1:t_1), wire_data(i,t_1)])/widths(i,wire);
            else
                nonlinear_energies(i,wire) = sum(squared_voltage...
                    - [wire_data(i,t_0), wire_data(i,t_0:t_1-1)].*[wire_data(i,t_0+1:t_1), wire_data(i,t_1)])/widths(i,wire);
            end
        else
            nonlinear_energies(i,wire) = 0;
        end
    end
end