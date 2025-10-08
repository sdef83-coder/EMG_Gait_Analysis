function n = get_cycle_count(META, pair)
    % Retourne un SCALAIRE = nombre de cycles sélectionnés pour la paire
    % Essaie de lire L_target_eff_<pair>, N_target_<pair>, target_<pair>
    v = try_get_meta_count(META, pair);
    n = 0;
    if isempty(v), return; end

    if islogical(v)
        n = sum(v(:));                 % masque logique -> somme
    elseif isnumeric(v)
        if isscalar(v)
            n = double(v);             % déjà un nombre
        else
            n = numel(v);              % vecteur / array d’indices -> longueur
        end
    elseif iscell(v)
        n = numel(v);                  % cell array -> nb d’éléments
    end
    if ~isfinite(n) || isnan(n)
        n = 0;
    end
end