function mm = molar_mass(el)
    % If el is like "Ca", will simply return molar mass of Ca.
    % If el is like ["Ca" 1; "O" 2] will return summed molar mass of CaO2.
    % data from https://www.angelo.edu/faculty/kboudrea/periodic/structure_mass.htm
    persistent els;
    persistent weights;
    if isempty(els)
        table = readtable('periodic_table_data.csv');
        els = string(table2array(table(:, 2)));
        weights = table2array(table(:, 4));
    end
    
    if length(el) == 1
        mm = weights(els==el);
    else
        counts = double(el(:, 2));
        els_ = el(:, 1);
        total = 0;
        for i=1:length(els_)
            partial_mm = weights(els==els_(i));
            total = total + counts(i) * partial_mm;
        end
        mm = total;
    end
end