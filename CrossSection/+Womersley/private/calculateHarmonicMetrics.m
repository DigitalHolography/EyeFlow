function h_metrics = calculateHarmonicMetrics(fitParams)
    arguments
        fitParams % list of fitParams based uppon the harmonics
    end

    fitParams = fitParams(:);
    leng_harmo = size(fitParams, 1);

    metrics = [fitParams.metrics];

    h_metrics = struct();

    % Shear harmonic ratios |τ2|/|τ1|, |τ3|/|τ1| (-)
    h_metrics.RhoTau21 = ternary_op(leng_harmo >= 2, @() norm(metrics(2).tau_n) / norm(metrics(1).tau_n), @() NaN);
    h_metrics.RhoTau31 = ternary_op(leng_harmo >= 3, @() norm(metrics(3).tau_n) / norm(metrics(1).tau_n), @() NaN);

    % Phase skewness of shear, ϕτ,2 − 2ϕτ,1 (°)
    h_metrics.DeltaPhiTau2 = ternary_op(leng_harmo >= 2, @() angle(metrics(2).tau_n) - 2 * angle(metrics(1).tau_n), @() NaN);

    % TODO: Add Description, because I can't Find anything
    h_metrics.RhoQ21       = ternary_op(leng_harmo >= 2, @() norm(metrics(2).Qn) / norm(metrics(1).Qn), @() NaN);
    h_metrics.RhoQ31       = ternary_op(leng_harmo >= 3, @() norm(metrics(3).Qn) / norm(metrics(1).Qn), @() NaN);
    h_metrics.DeltaPhiQ2   = ternary_op(leng_harmo >= 2, @() angle(metrics(2).Qn) - 2 * angle(metrics(1).Qn), @() NaN);
end