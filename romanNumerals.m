function [s] = romanNumerals(n)
%	By Francois Beauducel <beauducel@ipgp.fr>
%	  Institut de Physique du Globe de Paris

if ~isnumeric(n)
	error('N must be numeric array (scalar, vector or matrix).')
end

s = cell(size(n));

for k = 1:numel(n)
	m = max(floor((log10(n(k)) - log10(5000))/3) + 1,0);
	for i = m:-1:0
		ss = roman(fix(n(k)/10^(3*i)));
		if i == m
			s{k} = ss;
		else
			s{k} = ['(',s{k},')',ss];
		end
		n(k) = mod(n(k),10^(3*i));
	end
end

% converts to string if n is a scalar
if numel(n) == 1
	s = s{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=roman(n)
% this subfunction converts numbers up to 4999

r = reshape('IVXLCDM   ',2,5);	% the 3 last blank chars are to avoid error for n >= 1000
x = '';

m = floor(log10(n)) + 1;	% m is the number of digit

% n is processed sequentially for each digit
for i = m:-1:1
	ii = fix(n/10^(i-1));	% ii is the digit (0 to 9)
	% Roman numeral is a concatenation of r(1:2,i) and r(1,i+1)
	% combination with regular rules (exception for 4000 = MMMM)
	% Note: the expression uses REPMAT behavior which returns empty
	% string for N <= 0
	x = [x,repmat(r(1,i),1,ii*(ii < 4 | (ii==4 & i==4)) + (ii == 9) + (ii==4 & i < 4)), ...
		   repmat([r(2,i),repmat(r(1,i),1,ii-5)],1,(ii >= 4 & ii <= 8 & i ~= 4)), ...
		   repmat(r(1,i+1),1,(ii == 9))];
	n = n - ii*10^(i-1);	% substract the most significant digit
end



end