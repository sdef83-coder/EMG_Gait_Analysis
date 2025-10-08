function v = try_get_meta_count(META, pair)
    % Priorité : L_target_<pair> puis N_target_<pair> puis target_<pair>
    fields_try = { ...
        sprintf('L_target_eff_%s', pair), ...
        sprintf('N_target_%s', pair), ...
        sprintf('target_%s',   pair) ...
    };
    v = [];
    for k = 1:numel(fields_try)
        f = fields_try{k};
        if isfield(META, f)
            v = META.(f);
            return;
        end
    end
end