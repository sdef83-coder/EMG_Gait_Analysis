%% === COHERENCE AREA SIGNIFICATIVE (Aire sous la courbe) — Ajout dans DATA ===
% Charge Coherence_<Participant>.mat (script 1), calcule les aires significatives par condition,
% jambe, paire, sous-phase et bande (Alpha/Beta/Gamma), les ajoute dans DATA,
% et ré-enregistre le .mat.
% MODIFICATION: Applique le seuil de significativité pour calculer seulement l'aire significative

clear; close all; clc;

% ===== Paramètres =====
participant_id = 'CTL_63';  % <- adapte si besoin
coh_root = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(coh_root, 'ALL');

% Fichier à lire/écrire (même dossier)
infile  = fullfile(all_dir,  ['Coherence_' participant_id '.mat']);
outfile = infile;  % on ré-écrit par-dessus

% Bandes (Hz)
FreqBands = struct('Alpha',[8 12], 'Beta',[13 30], 'Gamma',[31 60]);

% ===== Chargement =====
if ~isfile(infile)
    error('Fichier introuvable: %s', infile);
end
L = load(infile);
if ~isfield(L, 'DATA')
    error('Le fichier ne contient pas la variable DATA: %s', infile);
end
DATA = L.DATA;

fprintf('>> Chargé: %s\n', infile);

% ===== Déduire les conditions valides =====
condNames = fieldnames(DATA);
keep = false(size(condNames));
for k = 1:numel(condNames)
    nm = condNames{k};
    keep(k) = isstruct(DATA.(nm)) && (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
end
condNames = condNames(keep);

sides = {'left','right'};
n_written = 0;
n_skipped = 0;
notes = {};

% ===== Boucles principales =====
for iC = 1:numel(condNames)
    condName = condNames{iC};
    fprintf('\n== Condition: %s ==\n', condName);

    % Meta informative
    try
        DATA.(condName).meta.CoherenceArea_FreqBands = FreqBands;
        DATA.(condName).meta.CoherenceArea_Method = ...
            'Mean over time (within phase), subtract significance threshold, set negative to 0, then trapz(Freq) over band. Freq sorted asc; NaNs filled.';
    end

    for s = 1:numel(sides)
        sideStr = sides{s};
        if ~isfield(DATA.(condName), sideStr), continue; end
        fprintf('-- Side: %s --\n', sideStr);
        Sside = DATA.(condName).(sideStr);

        % Fréquence requise
        if ~isfield(Sside, 'Freq') || isempty(Sside.Freq)
            notes{end+1} = sprintf('%s | %s : Freq manquant', condName, sideStr); 
            continue;
        end
        Freq = Sside.Freq(:);

        % Parcours des matrices "Coherence_*"
        fns = fieldnames(Sside);
        for iF = 1:numel(fns)
            fn = fns{iF};
            if ~startsWith(fn,'Coherence'), continue; end

            % Formats attendus:
            %  - Coherence_<m1>_<m2>
            %  - Coherence_<Phase>_<m1>_<m2>, Phase∈{LoadingResponse,MidStance,PreSwing,Swing}
            parts = strsplit(fn,'_');
            phase = ''; m1 = ''; m2 = '';
            if numel(parts) == 3
                % Coherence_<m1>_<m2>
                m1 = parts{2}; m2 = parts{3};
            elseif numel(parts) == 4 && any(strcmp(parts{2},{'LoadingResponse','MidStance','PreSwing','Swing'}))
                % Coherence_<Phase>_<m1>_<m2>
                phase = parts{2}; m1 = parts{3}; m2 = parts{4};
            else
                continue; % autre champ
            end

            % Chercher le seuil de significativité correspondant
            % Structure: DATA.condition.side.Seuil_muscle1_muscle2
            seuil_field = sprintf('Seuil_%s_%s', m1, m2);
            
            if ~isfield(Sside, seuil_field) || isempty(Sside.(seuil_field))
                notes{end+1} = sprintf('%s | %s | %s : Seuil manquant (%s)', ...
                    condName, sideStr, fn, seuil_field);
                n_skipped = n_skipped + 1;
                continue;
            end
            
            seuil = Sside.(seuil_field);

            C = Sside.(fn);
            if isempty(C) || size(C,1) ~= numel(Freq)
                notes{end+1} = sprintf('%s | %s | %s : dims mismatch (C=%s, nFreq=%d)', ...
                    condName, sideStr, fn, mat2str(size(C)), numel(Freq));
                continue;
            end

            % Spectre cohérence(f) = moyenne temporelle
            CohMean = mean(C, 2, 'omitnan');
            
            % ===== APPLICATION DU SEUIL DE SIGNIFICATIVITE =====
            % Soustraire le seuil à toutes les valeurs
            CohSignificant = CohMean - seuil;
            
            % Mettre à 0 toutes les valeurs négatives (non-significatives)
            CohSignificant(CohSignificant < 0) = 0;
            
            fprintf('   %s: Seuil=%.4f, Orig max=%.4f, Signif max=%.4f\n', ...
                fn, seuil, max(CohMean), max(CohSignificant));

            % Intégration par bande (sécuriser l'ordre + NaNs)
            [Fsort, idx] = sort(Freq, 'ascend');
            Csort = CohSignificant(idx);  % Utiliser la cohérence significative
            Csort = fillmissing(Csort,'linear');
            Csort = fillmissing(Csort,'nearest');

            bandNames = fieldnames(FreqBands);
            for ib = 1:numel(bandNames)
                bname  = bandNames{ib};
                brange = FreqBands.(bname);

                mask = (Fsort >= brange(1)) & (Fsort <= brange(2));
                area_val = NaN;
                if any(mask)
                    area_val = trapz(Fsort(mask), Csort(mask));
                end

                % Noms des champs pour aires significatives
                if isempty(phase)
                    outfield = sprintf('CoherenceArea_%s_%s_%s', bname, m1, m2);
                else
                    outfield = sprintf('CoherenceArea_%s_%s_%s_%s', phase, bname, m1, m2);
                end

                DATA.(condName).(sideStr).(outfield) = area_val;
                n_written = n_written + 1;
            end
        end
    end
end

fprintf('\n>> Aires écrites: %d champs.\n', n_written);
fprintf('>> Paires sautées (seuil manquant): %d\n', n_skipped);
if ~isempty(notes)
    fprintf('>> Notes (%d):\n', numel(notes));
    for k = 1:numel(notes), fprintf('   - %s\n', notes{k}); end
end

% ===== Sauvegarde =====
save(infile, 'DATA', '-v7.3');       % met à jour dans le dossier du participant
fprintf('>> Sauvegardé: %s\n', infile);