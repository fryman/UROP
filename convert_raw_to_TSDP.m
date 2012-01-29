function convert_raw_to_TSDP( file_without_extension )
%convert_raw_to_TSDP_FN Saves a .mat from a .ntt with TS and DP
%   Given a filename without its extension <file_without_extension>, we search for the .NTT in the
%   Data folder. Then we extract the TimeStamps and DataPoints from the
%   .NTT and save both variables as a .MAT file with the same name as the
%   <file_without_extension> arg passed.

fwe = file_without_extension;
pathToFoo = 'C:\Users\Alexander\Desktop\Aaron Fryman\Data'; %'F:\AaronF\Data\Manual\';
foo = fullfile(pathToFoo, strcat(fwe, '.ntt'));

[TimeStamps, DataPoints, Header] = Nlx2MatSpike_411(foo, [1 0 0 0 1], 1, 1, []);

line_adbits = Header{15};
adbit_chars = line_adbits(13:end-1);
adbit_vals = strread(adbit_chars);

number_of_samples = size(DataPoints,3);
points_in_waveform = 32;
scaled_dp_in_volts = DataPoints .* repmat(adbit_vals, [points_in_waveform 1 number_of_samples]);
scaled_dp_in_microvolts = scaled_dp_in_volts*10^6;
dp = reshape(scaled_dp_in_microvolts, 128, []);

savedFoo = fullfile(pathToFoo, fwe);
save(savedFoo, 'TimeStamps','dp');

end