function B = efficient_Ak_dotstar_A(A,k,strin,varargin)
% EFFICIENT_AK_DOTSTAR_A Compute A^k .* A in a memory efficient manner
%
% B = efficient_Ak_dotstar_A(A,k) computes B = (A^k).*A in a memory
% efficient manner by using a block of columns at a time.
%
% ... = efficient_Ak_dotstar_A(...,'key',value,'key',value) provides
% optional arguments that control the behvaior
%
%   'memory' : use up to this many bytes of memory (ignored if memory is
%     less than one dense vector of size A), the default is the size of A
%     itself in bytes, so the entire algorithm will take 3 times the memory
%     of the matrix A. This can be specified by "-1", listing 0 uses
%     minimal memory.
%   'disp' : display every k seconds, default is a doubling, which can be
%     specified by "-1", 0 turns off all displays

% Nicole Eikmeier and David F. Gleich
% Purdue University, 2015

assert(size(A,1) == size(A,2));

p = inputParser;
p.addOptional('memory',-1,@isnumeric);
p.addOptional('disp',-1,@isnumeric);
p.parse(varargin{:})
opts = p.Results;

n = size(A,1);
Adata = whos('A');

minmem = n*8;

if nargin < 3
strin = 'placeholder';
end

if opts.memory == -1
    blockmem = max(Adata.bytes,minmem);
else
    blockmem = max(minmem,opts.memory);
end

blocksize = floor(blockmem/(8*n)); % this is the number of columns
nblocks = ceil(n/blocksize);

% Now handle the display
nextstatus = opts.disp;
if nextstatus == 0
    nextstatus = Inf; % never report
elseif nextstatus < 0
    % use our geometric progession
    nextstatus = 15;
end    
t0 = tic;

Bcell = cell(1,nblocks);
for bi=1:nblocks
    blockregion = ((bi-1)*blocksize+1):min((bi)*blocksize,n);
    C = full(A(:,blockregion));
    for i=2:k
        C = A*C;
    end
    Bcell{bi} = sparse(C.*A(:,blockregion));
    if toc(t0) >= nextstatus
        % make a report!
        fprintf('Ak_dot_A : %9i of %9i columns done (%7.1f sec, %5.2f%% complete)\n', ...
            max(blockregion), n, toc(t0), 100*max(blockregion)/n);
        if opts.disp > 0
            nextstatus = nextstatus + opts.disp;
        else
            nextstatus = nextstatus + min(3600,2*nextstatus); % max of one hour
        end
    end
end
B = cell2mat(Bcell);

l = strfind(strin,'placeholder');
if isempty(l)
save(strin,'B');
end
%save(strin,'B');