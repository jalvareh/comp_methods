function spectrogram(w, type, plotmovie)

    close all;

    load handel
    v = y'/2;
    L = Fs;
    n = length(v);
    t2 = linspace(0,L,n+1);
    t = t2(1:n)/1000;
    k = (2*pi/L)*[0:n/2 -n/2:-1];
    ks = fftshift(k);
    kmax = max(ks);
    S = v;
    St = fft(S);
    
    if strcmp(type, 'piano');     
         figure(1);
         tr_piano=16;  % record time in seconds
         y=wavread('music1'); Fs=length(y)/tr_piano;
         plot((1:length(y))/Fs,y);
         xlabel('Time [sec]'); ylabel('Amplitude');
         title('Mary had a little lamb (piano)');  drawnow
         p8 = audioplayer(y,Fs); playblocking(p8);
    elseif strcmp(type, 'recorder')   
         figure(2);
         tr_rec=14;  % record time in seconds
         y=wavread('music2'); Fs=length(y)/tr_rec;
         plot((1:length(y))/Fs,y);
         xlabel('Time [sec]'); ylabel('Amplitude');
         title('Mary had a little lamb (recorder)');
         p8 = audioplayer(y,Fs); playblocking(p8);
    else
        return
    end

    % The following code makes a dynamical movie of this process as the 
    % parameter b is translated.
%     figure(1);
%     tslide=0:0.1:10;
%     spec = {};
%     for k = 1:length(w)
% 
%         Sgt_spec=[];
%         for j=1:length(tslide/2)
%             
%             % Gabor transform
%             tt = (t-tslide(j));
%             if strcmp(type,'gauss')
%                 g = exp(-w(k)*tt.^2); 
%             elseif strcmp(type,'haar')
%                 g = zeros(1,length(t));
%                 g(0 <= w(k)*tt & w(k)*tt < 1/2) = 1;
%                 g(1/2 <= w(k)*tt & w(k)*tt < 1) = -1;
%             elseif strcmp(type,'hat')
%                 g = (1-w(k)*tt.^2).*exp(-w(k)*tt.^2);
%             else
%                 break
%             end
%                         
%             
%             Sg=g.*S; Sgt=fft(Sg);
%             Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
%             
%             if plotmovie
%                 figure(1);
%                 subplot(3,1,1), plot(t,S,'k',t,g,'r'), ylabel('S(t)'), xlabel('Time (t)');
%                 axis([0 max(t) -1 1]);
%                 subplot(3,1,2), plot(t,Sg,'k'), ylabel('g*S(t)'), xlabel('Time (t)');
%                 axis([0 max(t) -kmax/100 kmax/100]);
%                 subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)));
%                 ylabel('FFT(g*S)'), xlabel('Frequency (\omega)');
%                 axis([0 kmax 0 1.1]);
%                 drawnow;
%             end
%             
%         end
% 
%         spec{k} = Sgt_spec;
%     end
%     
%     % The code just developed produces a matrix Sgt spec which contains the 
%     % Fourier transform at each slice in time of the parameter b. It is this
%     % matrix that produces the spectrogram of the time-frequency signal.
%     
%     if length(w) < 4
%         plotmax = length(w);
%     else
%         plotmax = 4;
%     end
%     
%     for k = 1:plotmax
% 
%         figure(2);
%         subplot(2,2,k);
%         pcolor(tslide,ks,spec{k}.'), shading interp, colormap(hot);
%         set(gca,'Ylim',[-kmax kmax],'Fontsize',14);
%         xlabel('Time (t)'); ylabel('Frequency (\omega)');
%         mystring = sprintf('Width = %0.1e s^{-1}',w(k));
%         title(mystring);
%         drawnow;
% 
%     end

    
end    




