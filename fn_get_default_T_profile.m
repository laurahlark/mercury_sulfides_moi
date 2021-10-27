function T_of_r=fn_get_default_T_profile()
    k2018_T_data = sortrows(load('knibbe_2018_ref_temp.csv'));
    T_of_r.r = k2018_T_data(:, 1)*1e5;
    T_of_r.T = k2018_T_data(:, 2);
end