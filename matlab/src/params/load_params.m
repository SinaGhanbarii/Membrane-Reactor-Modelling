function p = load_params()
%LOAD_PARAMS  Read data/params_literature.json (cached after the first call).
%   Every leaf is a struct with fields value, unit, source, kind.
%   A JSON null becomes [] here; use pget() to read a value, because pget
%   stops with an error when the value is missing.
persistent cache
if isempty(cache)
    here = fileparts(mfilename('fullpath'));          % .../matlab/src/params
    f = fullfile(here, '..', '..', 'data', 'params_literature.json');
    if ~isfile(f)
        error('load_params:missing', 'Data file not found: %s', f);
    end
    cache = jsondecode(fileread(f));
end
p = cache;
end
