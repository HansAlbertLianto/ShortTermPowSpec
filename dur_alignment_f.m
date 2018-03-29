
function [outputDataPath1, outputDataPath2] = dur_alignment_f(inputDataPath1, inputDataPath2)
%%%%%%%%%%%%%%%%%%%% Duration alignment code by PhD Student in MICL %%%%%%%%%%%%%%%%%%%%%%%

% Function call adds scripts to MATLAB Path so that the relevant scripts
% and functions can run.
full_path = mfilename('fullpath');
cropped_full_path = full_path(1:(length(full_path)-16));
tool_path = strcat(full_path(1:(length(full_path)-16)), '\lib');
addpath([cropped_full_path '\Processing']);
addpath([tool_path '\STRAIGHTV40_007d_512']);
addpath([tool_path '\DTW']);
addpath([tool_path '\fbank']);
addpath([tool_path '\spl_v3']);

% Default parameter initialization.
prm.defaultFrameLength = 4.0;
prm.F0frameUpdateInterval=10.000000;
prm.F0archUpperBound=600;
prm.F0archLowerBound=50;
prm.spectralUpdateInterval=10.000000;

% Frame size and frame overlaps.
f_size = 25;
f_overlap = 15;

% Resample all .wav files to 16000 samples per second.
fs_ref = 16000;

for i = length(inputDataPath1):-1:1
    index1 = i;
    if ((inputDataPath1(i) == '\') || (inputDataPath1(i) == '/'))
        break
    end
end
    
for j = length(inputDataPath2):-1:1
    index2 = j;
    if ((inputDataPath2(j) == '\') || (inputDataPath2(j) == '/'))
        break
    end
end  

% Getting the name of the file and the file extension.
source_file = inputDataPath1((index1 + 1):length(inputDataPath1));
target_file = inputDataPath2((index2 + 1):length(inputDataPath2));

% Initializing a features folder to align features.
feature_path = './features';
if ~exist(feature_path)
    mkdir(feature_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source_wav = inputDataPath1;
target_wav = inputDataPath2;

% Extract features using the STRAIGHT Library from both the source and
% target file.

source_feature = sprintf('%s/%s_ori.mat', feature_path, source_file);
if 1 %~exist(source_feature)
    [x,fs]=audioread(source_wav);
    if fs ~= fs_ref
        [p, q] = rat(fs_ref/fs);
        % disp(p);
        % disp(q);
        x = resample(x,p,q);
        fs = fs_ref;
    end
    [f0_s, ap_s] = exstraightsource(x,fs, prm);
    [sp_s] = exstraightspec(x, f0_s, fs, prm);
    
    frame_s = buffer(x, fs*f_size/1000, fs*f_overlap/1000);
    
    intensity_s = rms(frame_s);
    
    pointspl = 20*log10(intensity_s/20e-6);
    
    len_f0 = length(f0_s);
    len_int = length(intensity_s);
    
    if len_int > len_f0
        intensity_s = intensity_s(1:len_f0);
    elseif len_int < len_f0
        intensity_tmp = repmat(intensity_s(end), 1, len_f0-len_int);
        intensity_s = [intensity_s, intensity_tmp];
    end
    
    save(source_feature, 'f0_s', 'ap_s', 'sp_s', 'fs', 'intensity_s', '-mat');
else
    load(source_feature);
end


target_feature = sprintf('%s/%s_ori.mat', feature_path, target_file);
if 1 % ~exist(target_feature)
    [y,fs]=audioread(target_wav);
    if fs ~= fs_ref
        [p, q] = rat(fs_ref/fs);
        disp(p);
        disp(q);
        y = resample(y,p,q);
        fs = fs_ref;
    end
    [f0_t, ap_t] = exstraightsource(y, fs, prm);
    [sp_t] = exstraightspec(y, f0_t, fs, prm);
    
    frame_t = buffer(y, fs*f_size/1000, fs*f_overlap/1000);
    intensity_t = rms(frame_t).^2;
    
    len_f0 = length(f0_t);
    len_int = length(intensity_t);
    
    if len_int > len_f0
        intensity_t = intensity_t(1:len_f0);
    elseif len_int < len_f0
        intensity_tmp = repmat(intensity_t(end), 1, len_f0-len_int);
        intensity_t = [intensity_t, intensity_tmp];
    end
    
    save(target_feature, 'f0_t', 'ap_t', 'sp_t', 'fs', 'intensity_t', '-mat');
else
    load(target_feature);
end

fs = 16000;

% Perform Dynamic Time Warping on the mentioned features.

[wts, ~] = fft2melmx(1024, 16000, 50, 1, 133.33, 6855.5, 1);
wts_matrix = wts(:, 1:513);

fbank_s = wts_matrix * sp_s;
fbank_t = wts_matrix * sp_t;

[temp_p,temp_q,D] = DynamicTimeWarping(fbank_s, fbank_t);

% SM = simmx(sp_s,sp_t);
% [temp_p,temp_q,C] = dp(1-SM);

cost_vector = zeros(1,length(temp_p));
for k = 1:length(temp_p)
    cost_vector(k) = sum((fbank_s(:,temp_p(k)) - fbank_t(:,temp_q(k))).^2);
end

% Align the target .wav file to the source .wav file.
[~, IX] = sort(cost_vector);
new_temp_p = temp_p(IX);
new_temp_q = temp_q(IX);
[CC,IA,IC] = unique(new_temp_p', 'first');
new_temp_p = new_temp_p(IA);
new_temp_q = new_temp_q(IA);
[temp_p, IX] = sort(new_temp_p);
temp_q = new_temp_q(IX);

% Resynthesis so the new aligned waveforms are formed.
f0_s = f0_s(temp_p);
ap_s = ap_s(:, temp_p);
sp_s = sp_s(:, temp_p);
intensity_s = intensity_s(temp_p);
vuv_s = find(f0_s>0);

f0_t = f0_t(temp_q);
ap_t = ap_t(:, temp_q);
sp_t = sp_t(:, temp_q);
intensity_t = intensity_t(temp_q);
vuv_t = find(f0_t>0);

align_index = intersect(vuv_s, vuv_t);

align_f0_s = f0_s(align_index);
align_f0_t = f0_t(align_index);

align_intensity_s = intensity_s(align_index);
align_intensity_t = intensity_t(align_index);

% figure; plot(f0_t(align_index));
% hold on; plot(f0_s(align_index));
% 
% figure; plot(diff(f0_t(align_index)));
% hold on; plot(diff(f0_s(align_index)));

mkdir(strcat(inputDataPath1(1:index1), 'aligned'));
mkdir(strcat(inputDataPath2(1:index2), 'aligned'));

outputDataPath1 = strcat(inputDataPath1(1:index1), 'aligned\', inputDataPath1((index1 + 1):length(inputDataPath1)));
outputDataPath2 = strcat(inputDataPath2(1:index2), 'aligned\', inputDataPath2((index2 + 1):length(inputDataPath2)));

% Write aligned audio data into the .wav files.

sx = exstraightsynth(f0_s, sp_s, ap_s, fs, prm);
aligned_source_file = outputDataPath1;
audiowrite( aligned_source_file, sx/max(abs(sx)), 16000 );

sy = exstraightsynth(f0_t, sp_t, ap_t, fs, prm);
aligned_target_file = outputDataPath2;
audiowrite( aligned_target_file, sy/max(abs(sy)), 16000 );

% Save aligned features in the features folder.

source_feature_m = sprintf('%s/%s_aligned.mat', feature_path, source_file);
save(source_feature_m, 'f0_s', 'ap_s', 'sp_s', 'fs', 'align_intensity_s', '-mat');

target_feature_m = sprintf('%s/%s_aligned.mat', feature_path, target_file);
save(target_feature_m, 'f0_t', 'ap_t', 'sp_t', 'fs', 'align_intensity_t', '-mat');

% figure; 
% plot(sx/max(abs(sx))); hold on;
% plot(sy/max(abs(sy)), 'r');
end
