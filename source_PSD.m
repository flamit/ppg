function [source_signals,pow_peakiness] = source_PSD(source_signals,...
    pow_peakiness,data,start,colour_title,fs,LF,HF,f_ica,do_plot)
%compute source signals (using fastICA) - calculate PSD estimates
%rate the PSD peakiness of source signals

[num_sources,len] = size(data);

try
    %fast ICA
    if f_ica == 1
        crgb = fastica(data);        
        %plots
        if do_plot == 1
            figure
            for source_no = 1:num_sources
                subplot(num_sources,1,source_no); plot(crgb(source_no,:))
                if source_no == 1
                    title(sprintf('Source signals after applying ICA: %s',...
                        colour_title))
                end
            end
        end
    %no fast ICA    
    else
        crgb = data;
    end
    source_signals(start:start+num_sources-1,:) = crgb;
    
    %rate the PSD peakiness of (source) signals
    if do_plot == 1
        figure
    end
    for source_no = 1:num_sources
        %measure peakiness of each signal
        source = crgb(source_no,:);
        [peakiness,~,~,~,f,pxx] = PSD_peak(source,fs,LF,HF);
        pow_peakiness(start+source_no-1) = peakiness;
        %plot graph
        if do_plot == 1
            subplot(num_sources,1,source_no); plot(f,pxx)
            xlabel('Frequency (Hz)')
            ylabel('Power (dB)')
            graph_title = 'Power spectral density estimates of source signals: ';
            if source_no == 1
                title(sprintf('%s%s',graph_title,colour_title))
            end
        end
    end
    
catch
    pow_peakiness(start:start+num_sources-1) = -1;
    source_signals(start:start+num_sources-1,:) = zeros(num_sources,len);
    %close window if unsuccessful
    close
end