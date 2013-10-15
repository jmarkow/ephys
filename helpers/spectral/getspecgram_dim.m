function [T,F,MINF,MAXF]=getspecgram_dim(NSAMPLES,N,NOVERLAP,NFFT,FS,MINF,MAXF)
%
%
%
%
%
%

% min and max f of the F vector

if nargin<7 | isempty(MAXF), MAXF=10e3; end
if nargin<6 | isempty(MINF), MINF=1e3; end

% get the row and column indices

if mod(NFFT,2)==0
	rows=(NFFT/2+1);
else
	rows=(NFFT+1)/2;
end

% get the column indices, hop size is n-overlap

columns=fix((NSAMPLES-NOVERLAP)/(N-NOVERLAP));

% frequency identities are linearly spaced to Nyquist

F=((0:rows-1)./(rows-1)).*(FS/2);

colidx=1+(0:(columns-1))*(N-NOVERLAP);

% frequency spacing starts at N/2, then continues by hop size
% in samples, finally convert to time by /FS

T=((colidx-1)+((N/2)'))/FS;

% finally get the min and max frequency indices

MINF=max(find(F<=MINF));

if isempty(MINF)
	MINF=1;
end

MAXF=min(find(F>=MAXF));

if isempty(MAXF)
	MAXF=length(F);
end
