function p = load_params()
%LOAD_PARAMS  Read data/params_register.csv into a nested struct (cached after the first call).
%   One CSV row per parameter.  Columns: id, value, unit, source, kind, note, vmin, vmax,
%   related, sweep_lo, sweep_hi.
%   kind = reported | calculated | assumed | inferred | placeholder | sweep
%        -> leaf struct with fields value, unit, source, kind, note, related, sweep
%           value  : [] when blank (= not available); a;b;c becomes a numeric column vector
%           related: text "type:target;type:target" (reuses, conflicts_with, needs,
%                    inferred_from, range_from, maybe_affected_by)
%           sweep  : [lo; hi] when both bounds are given, otherwise []
%   kind = range -> numeric [vmin; vmax]   (validity range of a rate law, membrane or inhibition set)
%   kind = meta  -> plain text (labels, notes)
%   Read a number with pget(), which stops on a missing value.
persistent cache
if isempty(cache)
    here = fileparts(mfilename('fullpath'));          % .../matlab/src/params
    f = fullfile(here, '..', '..', 'data', 'params_register.csv');
    if ~isfile(f)
        error('load_params:missing', 'Data file not found: %s', f);
    end
    opts = detectImportOptions(f, 'Delimiter', ',', 'TextType', 'string');
    opts = setvartype(opts, opts.VariableNames, 'string');
    T = readtable(f, opts);
    need = ["id","value","unit","source","kind","note","vmin","vmax","related","sweep_lo","sweep_hi"];
    if ~all(ismember(need, string(T.Properties.VariableNames)))
        error('load_params:columns', 'params_register.csv must have columns: %s', strjoin(need, ', '));
    end
    if numel(unique(T.id)) ~= height(T)
        [~, i] = unique(T.id); dup = T.id(setdiff(1:height(T), i));
        error('load_params:duplicate', 'Duplicate parameter id: %s', strjoin(dup, ', '));
    end
    p = struct();
    for k = 1:height(T)
        parts = cellstr(split(T.id(k), '.'))';
        if ~all(cellfun(@isvarname, parts))
            error('load_params:badid', 'Row %d: "%s" is not a valid MATLAB field path.', k, T.id(k));
        end
        switch T.kind(k)
            case "range"
                lo = str2double(T.vmin(k)); hi = str2double(T.vmax(k));
                if isnan(lo) || isnan(hi)
                    error('load_params:badrange', 'Row %d (%s): vmin/vmax must be numbers.', k, T.id(k));
                end
                v = [lo; hi];
            case "meta"
                v = txt(T.value(k));
            otherwise
                sw = [];
                lo = str2double(T.sweep_lo(k)); hi = str2double(T.sweep_hi(k));
                if ~isnan(lo) && ~isnan(hi)
                    if lo >= hi, error('load_params:badsweep', 'Row %d (%s): sweep_lo must be below sweep_hi.', k, T.id(k)); end
                    sw = [lo; hi];
                end
                v = struct('value', parse_value(T.value(k), T.id(k)), 'unit', txt(T.unit(k)), ...
                           'source', txt(T.source(k)), 'kind', txt(T.kind(k)), 'note', txt(T.note(k)), ...
                           'related', txt(T.related(k)), 'sweep', sw);
        end
        p = setfield(p, parts{:}, v); %#ok<SFLD>
    end
    cache = p;
end
p = cache;
end

function s = txt(x)
if ismissing(x), s = ''; else, s = char(x); end
end

function v = parse_value(x, id)
if ismissing(x) || strlength(strtrim(x)) == 0, v = []; return; end
if contains(x, ';')
    v = str2double(split(x, ';'));
    if any(isnan(v)), error('load_params:badlist', '%s: list "%s" is not all numeric.', id, x); end
else
    n = str2double(x);
    if isnan(n), v = char(x); else, v = n; end
end
end
