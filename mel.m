function m=mel(f)
% converts to Mel scale

  m=2595.*log10(1+f./700);
