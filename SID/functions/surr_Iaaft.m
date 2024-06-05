%% iterative amplitude adjusted fourier tranform surrogates 
% algorithm of Schreiber and Schmitz - Physical Review Letters 1996

% y: series to replace
% nit: Number of Iteractions (default 7)
% stop: 'spe' preserves the spectrum, 'dis' preserves the distribution

function ys=surr_Iaaft(y,nit,stop)

narginchk(1,3);
if nargin < 3, stop='spe'; end %default to preserve the spectrum
if nargin < 2, nit=7; end %default 7 iterations

[ysorted,~]=sort(y); % from the lowest to the highest value
my=abs(fft(y));
ys=surr_ShuffColumn(y); % shuffling

for i=1:nit
    % step 1: impose the spectrum
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys); ys=real(ys);
    ys=ys-mean(ys);

    % step 2: impose the distribution
    [~,ysindice]=sort(ys);
    ypermuted=zeros(length(y),1);
    for ii=1:length(y)
        ypermuted(ysindice(ii))=ysorted(ii);
    end
    ys=ypermuted;

end

%imposing the same spectrum
if strcmp(stop,'spe')
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys);
end




