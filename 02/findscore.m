function findscore(w, type, plotmovie)

    close all;

    if strcmp(type, 'piano');     
         figure(1);
         tr_piano=16;  % record time in seconds
         y=audioread('music1.wav'); Fs=length(y)/tr_piano;
%          figure(1);
%          plot((1:length(y))/Fs,y);
%          xlabel('Time [sec]'); ylabel('Amplitude');
%          title('Mary had a little lamb (piano)');  drawnow
%          p8 = audioplayer(y,Fs); playblocking(p8);
    elseif strcmp(type, 'recorder')   
         figure(2);
         tr_rec=14;  % record time in seconds
         y=audioread('music2.wav'); Fs=length(y)/tr_rec;
%          figure(1);
%          plot((1:length(y))/Fs,y);
%          xlabel('Time [sec]'); ylabel('Amplitude');
%          title('Mary had a little lamb (recorder)');
%          p8 = audioplayer(y,Fs); playblocking(p8);
    else
        return
    end
    
    v = y'/2;
    L = Fs;
    n = length(v);
    t2 = linspace(0,L,n+1);
    t = t2(1:n)/1000;
    k = (2*pi/L)*[0:n/2-1 -n/2:-1];
    ks = fftshift(k);
    kmax = max(ks);
    S = v;
    St = fft(S);

%     The following code makes a dynamical movie of this process as the 
%     parameter b is translated.
    figure(1);
    tslide=0:0.1:10;
    spec = {};
    notes = {};
    for l = 1:length(w)

        Sgt_spec=[];
        notes{l} = zeros(1,length(tslide));
        
        for j=1:length(tslide)
            
            % Gabor transform
            tt = (t-tslide(j));
            g = exp(-w(l)*tt.^2);  
            
            Sg=g.*S; Sgt=fft(Sg);
            [dummy, idx] = max(abs(Sgt));
            note = 2*pi*abs(ks(idx))
            notes{l}(j) = note;
            Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
            
            if plotmovie
                figure(2);
                subplot(3,1,1), plot(t,S,'k',t,g,'r'), ylabel('S(t)'), xlabel('Time (t)');
                axis([0 max(t) -1 1]);
                subplot(3,1,2), plot(t,Sg,'k'), ylabel('g*S(t)'), xlabel('Time (t)');
                axis([0 max(t) -kmax/100 kmax/100]);
                subplot(3,1,3), plot(2*pi*ks,abs(fftshift(Sgt))/max(abs(Sgt)));
                ylabel('FFT(g*S)'), xlabel('\omega');
                axis([0 kmax 0 1.1]);
%                 subplot(4,1,4), plot(t(1:j),notes{l}(1:j),'k-');
%                 axis([0 max(t) 0 500]);
                drawnow;
            end
            
        end

        spec{l} = Sgt_spec;
    end
    
    % The code just developed produces a matrix Sgt spec which contains the 
    % Fourier transform at each slice in time of the parameter b. It is this
    % matrix that produces the spectrogram of the time-frequency signal.
    
    if length(w) < 4
        plotmax = length(w);
    else
        plotmax = 4;
    end
    
    for k = 1:plotmax

        figure(2);
        subplot(2,2,l);
        pcolor(tslide,2*pi*ks,spec{l}.'), shading interp, colormap(hot);
        set(gca,'Ylim',[0 500],'Fontsize',14);
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        mystring = sprintf('Width = %0.1e Hz',w(l));
        title(mystring);
        drawnow;

    end

    
end    




