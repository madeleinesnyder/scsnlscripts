function FilterDataMatrix = BandPassProc (DataMatrix)
TR = 2;
fl = 0.008;
fh = 0.1;
% Sampling Frequency
Fs = 1/TR;

% Set bandpass filter parameters
% Center Frequency
Fc = 0.5*(fh + fl);
% Filter Order
No = floor(Fs * 2/Fc);

disp('Filtering ........................................................');
[M,N,S,T] = size(DataMatrix);
FilterDataMatrix = zeros(size(DataMatrix));
% FIR filter Coefficients
B = getFiltCoeffs(zeros(1,T),Fs,fl,fh,0,No);
A = 1;

for s = 1:S
  % Extract time series of nth slice
  X = DataMatrix(:,:,s,:);
  Y = reshape(X,M*N,T)';
  % Filter the data
  Ya = filtfilt(B,A,Y);
  Yb = reshape(Ya',M,N,T);
  FilterDataMatrix(:,:,s,:) = Yb;
end
end