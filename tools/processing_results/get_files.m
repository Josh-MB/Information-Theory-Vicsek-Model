function [fname,E,R] = get_files(subdir)

fdir = fullfile(getenv('DATADIR'),'vicsek',subdir);

files = fullfile(fdir,'*.bin');
fprintf('finding files : %s\n\n',files);
fname = struct2cell(dir(files));
fname = fname(1,:)';
nfiles = length(fname);

E = zeros(nfiles,1);
R = zeros(nfiles,1);
for f = 1:nfiles
    sfname = fname{f};
    i(1:2) = strfind(sfname,'_');
    i(3) = strfind(sfname,'.bin');
    E(f) = str2double(sfname(i(1)+1:i(2)-1));
    R(f) = str2double(sfname(i(2)+1:i(3)-1));
    fname{f} = fullfile(fdir,fname{f});
end
assert(min(R) == 1);
numr = max(R);
nume = nfiles/numr;
assert(nume*numr == nfiles);
fname = reshape(fname,numr,nume)';

E = reshape(E,numr,nume);
R = reshape(R,numr,nume);
for r = 2:numr
    assert(isequal(E(r,:),E(1,:)));
end
for e = 2:nume
    assert(isequal(R(:,e),R(:,1)));
end
E = E(1,:);
R = R(:,1);
