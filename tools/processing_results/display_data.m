function I = display_data(fname,E,R,outputFile)

two_dim_hist = 0;

nume = length(E);
numr = length(R);

It = zeros(numr,nume);
H = zeros(numr,nume);
Iw = zeros(numr,nume);
Is = zeros(numr,nume);
Ilocal = zeros(numr,nume);
KLgte = zeros(numr,nume);
parm = read_data(fname,E,R,1,1);
U = 500;%parm.U
za = zeros(U,numr);
z = zeros(nume,1);
binders = zeros(numr,nume);
suscept = zeros(numr,nume);
order = zeros(numr,nume);
numReps = 7;
repBinWidth=numr/numReps;
if repBinWidth < 1
	repBinWidth = 1;
end
fullMI1DHist = zeros(parm.B,nume,numReps);
fullGTE1DHist = zeros(parm.B,nume,numReps);
fullGTE2DHist = zeros(parm.B,parm.B,nume,numReps);
fullMI1DHistMarginals = zeros(parm.B,nume);
zaccum = zeros(parm.U,repBinWidth, nume, numReps);
susceptGrouped = zeros(numReps, nume);
bindersGrouped = zeros(numReps,nume);

metric = 'Unknown';
metricNum = 0;
if isempty(strfind(fname{1,1},'mibin_')) == 0
	metric = 'Mutual Information';
	metricNum = 1;
elseif isempty(strfind(fname{1,1},'gtebin_')) == 0
	metric = 'Global Transfer Entropy';
	metricNum = 3;
elseif isempty(strfind(fname{1,1},'tebin_')) == 0
	metric = 'Transfer Entropy';
	metricNum = 2;
end

imethod = 'Unknown';
if parm.imethod == 0
	imethod = 'Metric';
elseif parm.imethod == 1
	imethod = 'Topological';
elseif parm.imethod == 2
	imethod = 'ISM';
end

umethod = 'Unknown';
if parm.umethod == 0
	umethod = 'Backward Update';
elseif parm.umethod == 1
	umethod = 'Forward Update';
end

for e = 1:nume
    for r = 1:numr
        [It(r,e),Iw(r,e),Is(r,e),Ilocal(r,e),H(r,e),ztest, terms, KLgte(r,e)] = calc_data(fname,E,R,e,r);
       
		rbin = floor((r-1)/repBinWidth)+1;
		% Accumulate ensemble statistics
		if metricNum == 1
			fullMI1DHist(:,e,rbin) = fullMI1DHist(:,e,rbin) + terms.count;
			fullMI1DHistMarginals(:,e) = fullMI1DHistMarginals(:,e) + terms.count2;
		end
		if metricNum == 3
			if parm.hist_gte_dims == 1
				fullGTE1DHist(:,e,rbin) = fullGTE1DHist(:,e,rbin) + terms.count;
			elseif parm.hist_gte_dims == 2
				fullGTE2DHist(:,:,e,rbin) = fullGTE2DHist(:,:,e,rbin) + terms.count;
			end
		end;
		order(r,e) = mean(ztest);

		idx = mod(r, repBinWidth)+1;
		zaccum(:,idx, e, rbin)=ztest(:);
		[suscept(r,e), binders(r,e)] = calcSusceptAndBinder(ztest);
    end
end

% Calculate ensemble order and susceptibility
for e = 1:nume
	for rbin = 1:numReps
		zacc = zaccum(:,:,e,rbin);
		[susceptGrouped(rbin, e), bindersGrouped(rbin, e)] = calcSusceptAndBinder(zacc);
	end
end

% Calculate ensemble statistics
MIHist = zeros(numReps, nume);
if metricNum == 1 && two_dim_hist == 0
	for r = 1:numReps
		for e = 1:nume
			P = fullMI1DHist(:,e,r)/sum(fullMI1DHist(:,e,r));
			% For discretised
			%MIHist(e) = (log(738.5826) - entropy(P))*log2(exp(1))
			%MIHist(e) = (log(4*parm.B/3) - (4/3)*entropy(P))*log2(exp(1))
			%MIHist(e) = (entropy(P))*log2(exp(1))
			% For continuous
			MIHist(r,e) = (log(parm.B) - entropy(P))*log2(exp(1));
		end
	end
end

GTEHist = zeros(numReps,nume);
if metricNum == 3 && two_dim_hist == 0
	if parm.hist_gte_dims == 1 
		for r = 1:numReps
			for e = 1:nume
				P = fullGTE1DHist(:,e,r)/sum(fullGTE1DHist(:,e,r));
				GTEHist(r,e) = (log(2*pi/(parm.B * E(e))) + entropy(P))*log2(exp(1));
			end
		end
	elseif parm.hist_gte_dims == 2
		for r = 1:numReps
			for e = 1:nume
				P = fullGTE2DHist(:,:,e,r)/sum(sum(fullGTE2DHist(:,:,e,r)));
				GTEHist(r,e) = (log(4*pi*pi/(parm.B*parm.B * E(e))) + centropy(P,2))*log2(exp(1));	
			end
		end
	end
end

% Plot figures
figure;
ax1=gca;
hold(ax1, 'on');

h = [];
if metricNum == 1 && any(MIHist(:))
	h(end+1) = plotErrorBar(ax1, E, MIHist, '-ro', '1D Histogram', numr);
elseif metricNum == 3 && any(GTEHist(:))
	if parm.hist_gte_dims == 1
		h(end+1) = plotErrorBar(ax1, E, GTEHist, '-ro', 'GTE 1D Histogram', numr);
	elseif parm.hist_gte_dims == 2
		h(end+1) = plotErrorBar(ax1, E, GTEHist, '-ro', 'GTE 2D Histogram', numr);
	end
elseif any(H)
	h(end+1) = plotErrorBar(ax1, E, H, '-bx', 'Histogram', numr);
end
h(end+1) = plotErrorBar(ax1, E, Iw, '-rx', 'NN', numr);
h(end+1) = plotErrorBar(ax1, E, Is, '-cx', 'NN - Decimation', numr);
h(end+1) = plotErrorBar(ax1, E, Ilocal, '-cx', 'NN - Localised', numr);
h(end+1) = plotErrorBar(ax1, E, KLgte, '-bo', '1D KL GTE', numr);
h(end+1) = plotErrorBar(ax1, E, order, '--kx', 'Order', numr);
h(end+1) = plotErrorBar(ax1, E, binders, '--mx', 'Binder Cumulant', numr);
h(end+1) = plotErrorBar(ax1, E, bindersGrouped, '--cx', 'Binder Cumulant (Grouped)', numr);

grid('on');

ax2 = axes('Position',get(ax1,'Position'),...
       'XAxisLocation','bottom',...
       'YAxisLocation','right',...
       'Color','none',...
       'XColor','k','YColor','k');
hold(ax2, 'on');
h(end+1) = plotErrorBar(ax2, E, suscept, '-go', 'Susceptibility', numr);
h(end+1) = plotErrorBar(ax2, E, susceptGrouped, '-cs', 'Susceptibility (Grouped)', numr);

h = nonzeros(h);

titleStr = sprintf('Vicsek (Cooling Regime) %s %s %s\nT_A=%d, N=%d, \\rho=%.2f, v=%.2f, T_k=%d, B=%d, %d reps\ngap=%d, \\chi=%.2f, J=%.2f, vis(\\eta)=%.2f, dt-factor=%.2f\0', imethod, umethod, metric, parm.U, parm.N, parm.rho, parm.v, parm.topo_neighbours, parm.B, numr, parm.record_T_steps,parm.chi,parm.J,parm.viscosity,parm.dt_factor);
title(titleStr);
xlabel('eta (noise)');
ylabel(ax1, metric);
ylabel(ax2, 'Susceptibility');
%legend('show', 'Location', 'southeast');
legend(h, 'Location', 'northeast');
%legend([h1, h2, h3], 'Location', 'northeast');
%legend([h1, h2, h3, h4, h5], 'Location', 'best');

hold off;

if nargin > 3
	pngFile = strcat(outputFile,'.png')
	figFile = strcat(outputFile,'.fig')
	print(gcf,pngFile,'-dpng','-r1000');
	figureHandle = gcf;
	saveas(figureHandle, figFile,'fig');
end
return
for e = 1:nume
	figure;
	bar(fullMI1DHistMarginals(:,e));
	fileName = sprintf('marginal_hist_eta_%.2f.png', E(e));
	print(gcf, fileName, '-dpng','-r1000');
end
end

function [handle] = plotErrorBar(ax, X, Y, lineSpec, name, numr)
	Y_m = mean(Y,1)';
	Y_s = std(Y,0,1)'/sqrt(numr);
	handle = 0;
	if any(Y_m)
		handle = errorbar(ax, X, Y_m, Y_s, lineSpec, 'DisplayName', name);
	end
end

function [chi, binder] = calcSusceptAndBinder(order)
	fourth = order.^4;
	squared = order.^2;
	binder = 1 - (mean(fourth) / (3*mean(squared)^2));
	chi = mean(squared) - mean(order)^2;
end