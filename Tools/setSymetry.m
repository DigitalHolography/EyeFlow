function v_profile = setSymetry(v_profile)
v_profile = (v_profile + v_profile(end:-1:1,:)) /2;
end