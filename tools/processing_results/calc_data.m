function [It, Iw, Is, Ilocal, Hist,z, terms, KLgte] = calc_data(fname,E,R,e,r)
    

[parms,z,count, count2] = read_data(fname,E,R,e,r);
It = parms.It*log2(exp(1));
Iw = parms.Iw*log2(exp(1));
Is = parms.Is*log2(exp(1));
Ilocal = parms.Ilocal*log2(exp(1));
KLgte = parms.KLgte*log2(exp(1));

if isempty(strfind(fname{1,1}, 'mibin_')) == 0
	P = count/sum(count(:));
	%Hist = 2*entropy(sum(P)) - entropy(P)
	Hist = entropy(sum(P,1)) + entropy(sum(P,2)) - entropy(P);
	terms.count = count;
	terms.count2 = count2;

elseif isempty(strfind(fname{1,1}, 'gtebin_')) == 0
	% Pretty sure this is wrong
	if parms.hist_gte_dims == 1
		P = count/sum(count(:));
		PP = sum(P,3);
		Hist = centropy(PP,1)-centropy(P,1);
	end
	terms.count = count;

elseif isempty(strfind(fname{1,1}, 'tebin_')) == 0
	P = count/sum(count(:));
	PP = sum(P,3);
	Hist = centropy(PP,1)-centropy(P,1);
	terms.count = count;
else
	Hist = 0;
end

Hist = Hist*log2(exp(1));

