function [parms,z,count, count2] = read_data(fname,E,R,e,r)

assert(e >= 1 && e <= length(E));
assert(r >= 1 && r <= length(R));

count2=0;
newname = strrep(fname{e,r}, '.bin', '.log.I');
fprintf('Reading file %s : eta = %f, RUN = %d\n',newname,E(e),R(r));
fprintf('Reading file %s : eta = %f, RUN = %d\n',fname{e,r},E(e),R(r));
fid = fopen(newname);
fidBin = fopen(fname{e,r});
fseek(fidBin,0,'eof');
nbytes=ftell(fidBin);
fseek(fidBin,0,'bof');
[vals, c] = fscanf(fid,'%f');
assert(c <= 5 && c >= 0);
parms.H = 0;%vals(1);
parms.It = 0;
parms.Iw = 0.0;
parms.Is = 0.0;
parms.Ilocal = 0.0;
parms.KLgte = 0.0;
if c >= 1
	parms.It = vals(1); % Per timestep
end
if c >= 2 % Whole series
	parms.Iw = vals(2);%vals(2);
end
if c >= 3 % Subsets
	parms.Is = vals(3);
end
if c >= 4 % Local flock
	parms.Ilocal = vals(4);
end
if c >= 5 % 1D KL GTE
	parms.KLgte = vals(5);
end

fclose(fid);

assert(check_magic(fidBin));

[parms.N,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.rho, c] = fread(fidBin,1,'double'); assert(c == 1);
[parms.v,   c] = fread(fidBin,1,'double'); assert(c == 1);
[parms.eta, c] = fread(fidBin,1,'double'); assert(c == 1);
[parms.L,   c] = fread(fidBin,1,'double'); assert(c == 1);
[parms.B,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.U,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.ksg_gte_dims, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.hist_gte_dims, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.GK,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.GKavg,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.topo_neighbours,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.ksg_neighbours,   c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.gen_seed, c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.sim_seed, c] = fread(fidBin,1,'uint64'); assert(c == 1);
[parms.imethod, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.umethod, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.te_shuffle_dim, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.discretise, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.ksg_local, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.record_T_steps, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.chi, c] = fread(fidBin, 1, 'double'); assert(c == 1);
[parms.J, c] = fread(fidBin, 1, 'double'); assert(c == 1);
[parms.viscosity, c] = fread(fidBin, 1, 'double'); assert(c == 1);
[parms.dt_factor, c] = fread(fidBin, 1, 'double'); assert(c == 1);
[parms.force_dt, c] = fread(fidBin, 1, 'double'); assert(c == 1);
[parms.rotate_frame, c] = fread(fidBin, 1, 'int32'); assert(c == 1);
[parms.shuffle, c] = fread(fidBin, 1, 'int32'); assert(c == 1);

U = parms.U/parms.record_T_steps;
[z1,c] = fread(fidBin,2*U,'double'); assert(c == 2*U);
z1 = reshape(z1,2,U);
z = sqrt(sum(abs(z1).^2,1));
atEnd=ftell(fidBin);
% fprintf('Consumed %d of %d bytes! (Not done yet though)', atEnd, nbytes);

B=parms.B;
if isempty(strfind(fname{1,1}, 'mibin_')) == 0
	[count,c] = fread(fidBin,B,'uint64'); assert(c==B);
	[count2, c] = fread(fidBin, B, 'uint64'); assert(c==B);
% Do gte first since it would also trigger tebin code
elseif isempty(strfind(fname{1,1}, 'gtebin_')) == 0
	if parms.hist_gte_dims == 1
		[count,c] = fread(fidBin,B,'uint64'); assert(c==B);
	elseif parms.hist_gte_dims == 2
		[count,c] = fread(fidBin,B*B,'uint64'); assert(c==B*B);
		count = reshape(count, [B B]);
	end
elseif isempty(strfind(fname{1,1}, 'tebin_')) == 0
	[count,c] = fread(fidBin,B*B*B,'uint64'); assert(c==B*B*B);
	count = reshape(count, [B B B]);
elseif isempty(strfind(fname{1,1}, 'params_')) == 0
	count = [];
end

atEnd=ftell(fidBin);
if atEnd ~= nbytes
	fprintf('Data remaining in file!!! Only consumed %d of %d bytes!\n', atEnd, nbytes);
else
	fprintf('Consumed whole file!\n');
end

fclose(fidBin);

end
