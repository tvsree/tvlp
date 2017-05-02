function f=imel(m)
% converts from Mel scale
  
  f=(10.^(m./2595)-1).*700;
