function v_profile = setBoundariesZeros(v_profile)

v_profile = v_profile - (v_profile(1,:)+v_profile(end,:))/2;

end