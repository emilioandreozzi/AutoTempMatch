function [Heartbeats,tmpl,wtmpl,NCC,params] = AutoTempMatch(sig,fs,varargin)

% AUTOTEMPMATCH implements an automated template matching method for
% ECG-free heartbeat detection in cardiomechanical signals, such as 
% Seismocardiograms, Gyrocardiograms, Forcecardiograms, 
% Phonocardiograms (heart sounds) and pulse wave signals (e.g., acquired
% via PPG sensors or force/pressure sensors placed onto
% peripheral blood vessels).
%
% The function only requires the input signal and its sampling frequency,
% and provides the temporal locations of heartbeats detected within the
% signal.
%
%
% When using this function, please cite the following article:
%
% 1) Parlato, S., Centracchio, J., Esposito, D., Bifulco, P., & Andreozzi,
%    E. (2025). Fully Automated Template Matching Method for ECG-Free
%    Heartbeat Detection in Cardiomechanical Signals of Healthy and
%    Pathological Subjects. Physical and Engineering Sciences in Medicine.
%
%
% Additional details on the method and its performance are reported in the
% following articles:
%
% 2) Centracchio, J., Parlato, S., Esposito, D., Bifulco, P., & Andreozzi,
%    E. (2023). ECG-Free Heartbeat Detection in Seismocardiography Signals
%    via Template Matching. Sensors (Basel, Switzerland), 23(10), 4684.
%    https://doi.org/10.3390/s23104684   
% 
% 3) Parlato, S., Centracchio, J., Esposito, D., Bifulco, P., & Andreozzi,
%    E. (2023). Heartbeat Detection in Gyrocardiography Signals without
%    Concurrent ECG Tracings. Sensors (Basel, Switzerland), 23(13), 6200.
%    https://doi.org/10.3390/s23136200   
% 
% 4) Parlato, S., Centracchio, J., Esposito, D., Bifulco, P., & Andreozzi,
%    E. (2023). ECG-Free Heartbeat Detection in Seismocardiography and
%    Gyrocardiography Signals Provides Acceptable Heart Rate Variability
%    Indices in Healthy and Pathological Subjects. Sensors (Basel,
%    Switzerland), 23(19), 8114. https://doi.org/10.3390/s23198114    
% 
% 5) Centracchio, J., Parlato, S., Esposito, D., & Andreozzi, E. (2024).
%    Accurate Localization of First and Second Heart Sounds via Template
%    Matching in Forcecardiography Signals. Sensors (Basel, Switzerland),
%    24(5), 1525. https://doi.org/10.3390/s24051525   
%
%
%
%
% INPUTS:
% - sig:        input signal
% - fs:         sampling frequency of the input signal
%
% PARAMETERS:
% The values of some parameters of the algorithm can be modified by the
% user, which has to provide them as additional input arguments to the
% function in the following format: 'ParameterName',ParameterValue
%
% Parameter names and descriptions:
% - 'TimeWin':    time window (in seconds) used to select the template (default is 10 s) 
% - 'pre':        value of the median IBI for the interval preceding the envelope peak (default is 0.2 s)
% - 'post':       value of the median IBI for the interval preceding the envelope peak (default is 0.5 s)
% - 'envopts':    cell array of parameters for the envelope extraction (default is rectification + LP filtering at 3 Hz, with MinPeakProminence of 0.25)
%                 Format of cell array:  {envelope function handle, MinPeakProminence for envelope peaks detection}
% - 'NCCpkprom':  minimum peak prominence for NCC peaks localization (default is 0.5)
% - 'NCCpkdist':  minimum peak distance for NCC peaks localization (default is 0.5 s)
% 
% 
% OUTPUTS:
% - Heartbeats: array of indices of localized heartbeats
% - tmpl:       selected template
% - wtmpl:      indices of samples of the selected template, referred to the whole input signal
% - NCC:        normalized cross-correlation function computed between the input signal "sig" and the selected template "tmpl"
% - params:     structure with actual parameters used for the algorithm
% 
% 
% SINTAX EXAMPLES:
%
% - BASIC:
%     Heartbeats = AutoTempMatch(sig,fs);
%  
% - ADVANCED:
%     Parameters: template selection via Hilbert envelope function, with
%     minimum peak prominence of 0.4 for envelope peaks detection, on 
%     5-second segments, with pre = 0.15 s and post = 0.70 s; NCC peaks
%     localization with minimum peak prominence of 0.7 for NCC peaks
%     localization.
%
%     envfunc = @(x) abs(hilbert(x));
%     [Heartbeats,tmpl,wtmpl,NCC,params] = AutoTempMatch(sig,fs,'TimeWin',5,'pre',0.15,'post',0.7,'envopts',{envfunc,0.4},'NCCpkprom',0.7);




%% PARAMETERS SETTING

% EXTRACTION OF PARAMETERS SPECIFIED BY THE USER

% Number of arguments
Narg = nargin-2;  % -2 to remove "sig" and "fs" which are not optional parameters and are considered in "nargin"

% Number of arguments must always be even, because the name and value of
% each parameter must be provided. So, if the number of arguments is odd,
% the user made an error, which must be reported.
if(mod(Narg,2) > 0)
    errstr = ['ERROR: wrong parameters passed as input arguments.\n' ...
    'Each parameter must be provided as a couple of ''ParameterName'',ParameterValue.' ...
    'Example: \n[Heartbeats] = AutoTempMatch(sig,fs,''TimeWin'',5,''NCCpkprom'',0.3);'];
    error(errstr);
end

% Extract the actual parameters specified by the user
for i=1:2:(Narg-1)
    eval(sprintf('%s = varargin{%d};',varargin{i},i+1))
end

clear i Narg varargin


% SETTING OF DEFAULT VALUES FOR PARAMETERS NOT SPECIFIED BY THE USER

% Time window for template selection
if(~exist('TimeWin','var'))
    TimeWin = 10;
else
    if(isempty(TimeWin))
        TimeWin = 10;
    end
end

% Pre and Post intervals for template selection
if(~exist('pre','var'))
    pre = 0.2;
else
    if(isempty(pre))
        pre = 0.2;
    end
end

if(~exist('post','var'))
    post = 0.5;
else
    if(isempty(post))
        post = 0.5;
    end
end

% Envelope for template selection
[b,a] = butter(2,3/(fs/2));
ENV = @(x) filtfilt(b,a,x.^4);
prom = 0.25;
if(~exist('envopts','var'))
    envopts = {ENV,prom};
else
    if(isempty(envopts))
        envopts = {ENV,prom};
    end   
end

% Minimum peak prominence for NCC peaks localization
if(~exist('NCCpkprom','var'))
    NCCpkprom = 0.5;
else
    if(isempty(NCCpkprom))
        NCCpkprom = 0.5;
    end   
end

% Minimum peak distance for NCC peaks localization
if(~exist('NCCpkdist','var'))
    NCCpkdist = fs/2;
else
    if(isempty(NCCpkdist))
        NCCpkdist = fs/2;
    end   
end

params.TimeWin = TimeWin;
params.pre = pre;
params.post = post;
params.envopts = envopts;
params.NCCpkprom = NCCpkprom;
params.NCCpkdist = NCCpkdist;


%% AUTOMATIC TEMPLATE SELECTION
% tmpl contains the selected template
% wtmpl contains the indices of the samples included in the selected
%       template, referred to the whole input signal "sig"
[tmpl,wtmpl] = AutoTemplSel(sig,fs,TimeWin,pre,post,envopts);

if(isempty(tmpl))
    fprintf('\n\nERROR:\nNo reliable segment was found for the selection of a heartbeat template.\n\n')
    Heartbeats = [];
    NCC = [];
else

    % Location of template absolute maximum
    [~,imax] = max(tmpl(1:(length(tmpl)/2)));
    
    
    
    %% NORMALIZED CROSS-CORRELATION
    NCC = nxcorr(sig,tmpl,imax); 
    NCC(NCC<0) = 0;
    cut = length(NCC)-length(sig);
    NCC = NCC(cut/2+1:end-cut/2); 
    
    
    
    %% NCC PEAKS LOCALIZATION
    [~,NCClocs] = findpeaks(NCC,'MinPeakProminence',NCCpkprom,'MinPeakDistance',NCCpkdist);
    Heartbeats = NCClocs;
end


end



function [tmpl,wtmpl] = AutoTemplSel(sig,fs,TimeWin,pre,post,envopts)
% This function performs an automatic selection of a heartbeat template
% from a cardio-mechanical signals, such as Seismocardiogram or
% Gyrocardiogram.

% The function first considers a time window at the beginning of the input
% signal, extracts the envelope of signal, locates the peaks 
% corresponding to heartbeats, and verify if at least 3 heartbeats have
% been identified in the time window, else it considers the next time
% window in the signal.

% If at least 3 heartbeats have been identified, the inter-beat intervals
% (IBI) are computed, and the median value is used to segment potential
% heartbeat templates.

% The starting point of a segment is obtained by considering a time
% interval preceding the envelope peak. The default value is 0.2, but the user
% can provide a different value as an argument of the function.

% The ending point of a segment is obtained by considering a time
% interval following the envelope peak. The default value is 0.5, but the user
% can provide a different value as an argument of the function.

% The function then computes the correlation coefficients between each
% couple of heartbeat segments and selects the one with the highest mean
% correlation with all others as the template.


% INPUTS:
% - sig: input signal
% - fs: sampling frequency of the input signal
% - TimeWin: time window (in seconds) used to select the template (default is 10 s) 
% - pre: value of the median IBI for the interval preceding the envelope peak
% - post: value of the median IBI for the interval preceding the envelope peak
% - envopts: array of parameters for the envelope extraction:
%            [envelope function handle, MinPeakProminence for envelope peaks detection]

% OUTPUTS:
% - tmpl: selected template
% - wtmpl: indices of samples of the selected template in the input signal

% SINTAX:
% - Basic:
%     tmpl = AutoTemplSel3(sig,fs);
%     10 seconds, pre: 0.20 s, post: 0.50 s, ENV: sig.^4 + 2nd ord Low-Pass
%     at 3 Hz cut-off, 0.25 MinPeakProminence for envelope peaks detection
%  
% - Advanced:
%     envfunc = @(x) abs(hilbert(x));
%     [tmpl,wtmpl] = AutoTemplSel3(sig,fs,5,0.15,0.7,{envfunc,0.4})
%     5 seconds, pre: 0.15 s, post: 0.70 s, ENV: hilbert envelope, 0.4
%     MinPeakProminence for envelope peaks detection


% %% DEFAULT OPTIONS
% 
% % TIME WINDOW
% if(~exist('TimeWin','var'))
%     TimeWin = 10;
% end
% 
% % PRE AND POST
% if(~exist('pre','var'))
%     pre = 0.2;
% end
% 
% if(~exist('post','var'))
%     post = 0.5;
% end
% 
% ENVELOPE
if(exist('envopts','var'))
    ENV = envopts{1};
    prom = envopts{2};
else
    [b,a] = butter(2,3/(fs/2));
    ENV = @(x) filtfilt(b,a,x.^4);
    prom = 0.25;
end



%% TEMPLATE SELECTION
normabs = @(x) x/max(abs(x(:)));
N = TimeWin*fs;
int2 = 1:N;

go = 0;
chunkOK = 0;

while(go == 0)
    
    % Current signal segment
    sigc = sig(int2);

    % ENVELOPE EXTRACTION  
    sigENV = ENV(sig);
    cut = sigENV(int2);

    % COARSE HEARTBEATS LOCALITAZION IN SIGNAL ENVELOPE
    [~,locsmax] = findpeaks(normabs(cut),'MinPeakProminence',prom);



    % RELIABILITY EVALUATION OF CURRENT SIGNAL CHUNK
    % The reliability is evaluated by assessing that:
    % a) a reasonable number of heartbeats have been recognized in the current chunk
    % b) the MAD of the inter-beat intervals of those heartbeats is below areasonable threshold

    % Minimum number of heartbeats is determined by considering a minimum
    % heart rate of 30 BPM 
    minbeats = floor(TimeWin/2);
    

    % SDNN threshold
    SDNN_thres = 0.335*fs;
    
    IBI = diff(locsmax);
    MAD = @(x) median(x-mean(x));
    MADNN = MAD(IBI);
    
    
    % Assessment of reliability criteria
    if( (numel(locsmax)>=minbeats) && (MADNN < SDNN_thres) )
        go = 1;
        chunkOK = 1;
    else
        int2 = int2(end)+1:(int2(end)+1+N);
        if(int2(end) > length(sig))
            go = 1;
            chunkOK = 0;
        else
            go = 0;
        end
    end

end

if(chunkOK == 1)
    
    % Heartbeats segmentation
    Npre = round(pre*fs);
    Npost = round(post*fs);
    
    clear w segs2
    k = 0;
    for i=1:numel(locsmax)
        if((locsmax(i)-Npre+1 >= 1) && locsmax(i)+Npost <= length(sigc))
            k = k+1;
            w(:,k) = locsmax(i)-Npre+1:locsmax(i)+Npost;
            segs2(:,k) = sigc(w(:,k));
        end
    end
    
    
    % Correlation coefficients matrix
    CM2 = maxXCmat(segs2);
    
    % Template selected as the heartbeat with the highest of mean/SD ratio of
    % correlation with other heartbeats 
    [~,iseg] = max(mean(CM2)./std(CM2));
    tmpl = segs2(:,iseg);
    wtmpl = w(:,iseg);
    
    % Adjust template window indices to refer them to the absolute time axis of the whole signal
    wtmpl = wtmpl + int2(1);
else
    tmpl = [];
    wtmpl = [];

end

end

function CM = maxXCmat(segs)

for i=1:size(segs,2)
    for j=1:size(segs,2)
        if(i==j)
            CM(i,j) = 1;
        else
            CM(i,j) = max(xcorr(segs(:,i),segs(:,j),'coeff'));
        end
    end
end

end