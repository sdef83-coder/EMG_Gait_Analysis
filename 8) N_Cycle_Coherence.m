%% === L_target_eff PAR CONDITION (Participant x Condition x Pair) ===
% - Parcourt les fichiers Coherence_CTL_XX.mat
% - Paires: TAprox_TAdist, GM_SOL, RF_ST, VL_RF, GMED_VL, GMED_RF
% - Lit DATA.(cond).meta.L_target_eff_<pair>
% - Déduit le nombre de cycles (scalaire, vecteur, struct courante, etc.)
% - Ecrit N_Cycle_Coherence.csv dans le dossier Results\Coherence

clc; clear; close all;

IN_DIR  = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
OUT_DIR = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
OUT_CSV = fullfile(OUT_DIR, 'N_Cycle_Coherence.csv');

if ~exist(IN_DIR,'dir'), error('Dossier introuvable: %s', IN_DIR); end
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

extract_id = @(fn) string(regexprep(fn, '^Coherence_([^\.]+)\.mat$', '$1'));

T = table('Size',[0 4], ...
          'VariableTypes', {'string','string','string','double'}, ...
          'VariableNames', {'Participant','Condition','Pair','L_target_eff'});

files = dir(fullfile(IN_DIR, 'Coherence_*.mat'));
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    S = load(fpath, 'DATA'); if ~isfield(S,'DATA'), warning('DATA absent: %s', files(i).name); continue; end
    DATA = S.DATA; pid = extract_id(files(i).name);

    condNames = fieldnames(DATA);
    keep = false(size(condNames));
    for k=1:numel(condNames)
        nm = condNames{k};
        keep(k) = isstruct(DATA.(nm)) && (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right') || isfield(DATA.(nm),'meta'));
    end
    condNames = condNames(keep);

    for ic = 1:numel(condNames)
        cond = condNames{ic};
        meta_here = []; if isfield(DATA.(cond),'meta') && isstruct(DATA.(cond).meta), meta_here = DATA.(cond).meta; end

        % 1) lister les paires via meta: L_target_eff_<m1>_<m2>
        pair_list = {};
        if ~isempty(meta_here)
            mf = fieldnames(meta_here);
            mask = startsWith(mf, 'L_target_eff_');
            keys = mf(mask);
            for j=1:numel(keys)
                % "L_target_eff_<m1>_<m2>"
                pair_list{end+1} = extractAfter(keys{j}, 'L_target_eff_'); %#ok<AGROW>
            end
            pair_list = unique(pair_list);
        end

        % 2) fallback: déduire via L_<m1>_<m2> sur left/right
        if isempty(pair_list)
            tmp = {};
            for side = {'left','right'}
                s = side{1};
                if isfield(DATA.(cond), s)
                    fns = fieldnames(DATA.(cond).(s));
                    lm = startsWith(fns, 'L_');
                    tmp = [tmp, extractAfter(fns(lm), 'L_')']; %#ok<AGROW>
                end
            end
            pair_list = unique(tmp);
        end
        if isempty(pair_list), continue; end

        % 3) lire L_target_eff_<pair> (sinon L_left+L_right)
        for ip = 1:numel(pair_list)
            pair = pair_list{ip};
            val = NaN;
            fmeta = ['L_target_eff_' pair];
            if ~isempty(meta_here) && isfield(meta_here, fmeta)
                val = double(meta_here.(fmeta));
            else
                fL = ['L_' pair];
                vL = NaN; vR = NaN;
                if isfield(DATA.(cond),'left')  && isfield(DATA.(cond).left,  fL), vL = double(DATA.(cond).left.(fL));  end
                if isfield(DATA.(cond),'right') && isfield(DATA.(cond).right, fL), vR = double(DATA.(cond).right.(fL)); end
                if ~isnan(vL) || ~isnan(vR), val = nansum([vL, vR]); end
            end
            if ~isnan(val)
                T = [T; {pid, string(cond), string(pair), val}]; %#ok<AGROW>
            end
        end
    end
end

T = sortrows(T, {'Participant','Condition','Pair'});
writetable(T, OUT_CSV);
fprintf('OK: %s\n', OUT_CSV);

%% ===== Fonctions locales =====
function id = local_extract_id(fname)
% Extrait "CTL_XX..." depuis "Coherence_CTL_XX....mat"
    id = "";
    tok = regexp(fname, 'Coherence_(?<id>CTL_[^\.]+)\.mat$', 'names');
    if ~isempty(tok) && isfield(tok,'id')
        id = string(tok.id);
    end
end

function n = local_get_n_cycles(x)
% Déduit le nombre de cycles à partir de formats courants:
% - scalaire: pris tel quel
% - vecteur/logique: longueur ou nnz
% - struct: essaie des champs usuels, sinon champ unique
% - cell: longueur
    n = NaN;

    if isnumeric(x)
        if isscalar(x), n = double(x); else, n = numel(x); end
        return
    end
    if islogical(x)
        % si logique vecteur -> nombre de true
        if isscalar(x), n = double(x~=0); else, n = nnz(x); end
        return
    end
    if iscell(x)
        n = numel(x); return
    end
    if isstruct(x)
        candidates = {'Nb_Cycles','N','n','count','Valid_Cycles','valid_cycles', ...
                      'indices','idx','cycles','cycle_idx','valid_idx'};
        for c = 1:numel(candidates)
            f = candidates{c};
            if isfield(x, f)
                val = x.(f);
                if isnumeric(val) || islogical(val)
                    if isscalar(val), n = double(val); else, n = numel(val); end
                elseif iscell(val)
                    n = numel(val);
                end
                if ~isnan(n), return, end
            end
        end
        % dernier recours: struct à champ unique contenant un vecteur/scalaire
        fns = fieldnames(x);
        if numel(fns)==1
            val = x.(fns{1});
            if isnumeric(val) || islogical(val)
                if isscalar(val), n = double(val); else, n = numel(val); end
                return
            elseif iscell(val)
                n = numel(val); return
            end
        end
    end
end