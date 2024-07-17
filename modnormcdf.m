function [ y ] = modnormcdf( x )

y = normcdf(x);
y(y>.5)=1;


end

