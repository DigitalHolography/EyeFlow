function ShowOutputs(paths, output_dir)
% This function show multiple outputs from the foldermanagement drawerlist

N = length(paths);

for path_idx = 1:N
    split_path = strsplit(paths{path_idx}, '\');
    main_foldername = split_path{end};
    folder_name = strcat(main_foldername, '_EF');
    ef_path = fullfile(paths{path_idx}, 'eyeflow');
    list_dir = dir(ef_path);
    idx = 0;

    for i = 1:length(list_dir)

        if contains(list_dir(i).name, folder_name)
            match = regexp(list_dir(i).name, '\d+$', 'match');

            if ~isempty(match) && str2double(match{1}) >= idx
                idx = str2double(match{1}); %suffix
            end

        end

    end

    last_folder_name = sprintf('%s_%d', folder_name, idx);

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'mask', [main_foldername, '_vesselMap.png']))
        segmentation_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'mask', [main_foldername, '_vesselMap.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_allrad_Artery_time.png']))
        bvr_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_allrad_Artery_time.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_f_artery_graph.png']))
        Arteries_fRMS_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_f_artery_graph.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_RI_velocityArtery.png']))
        ARI_velocity_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_RI_velocityArtery.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_histogramVelocityArtery.png']))
        histo_art_velocity_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_histogramVelocityArtery.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_strokeAndTotalVolume_Artery.png']))
        Stroke_total_volume{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_strokeAndTotalVolume_Artery.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_vessels_velocity_graph.png']))
        Vessels_velocity{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_vessels_velocity_graph.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_RI_velocityVein.png']))
        VRI_velocity_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_RI_velocityVein.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_A_sections.png']))
        A_sections{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_A_sections.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_diasys_Artery.png']))
        diasys_Artery{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_diasys_Artery.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_diasys_Vein.png']))
        diasys_Vein{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'crossSectionsAnalysis', [main_foldername, '_diasys_Vein.png']);
    end

    if isfile(fullfile(ef_path,d last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_ArterialWaveformAnalysis_v_artery.png']))
        ArterialWaveformAnalysis_artery{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_ArterialWaveformAnalysis_v_artery.png']);
    end

end

[l, L] = bestMontageLayout(N);

figure(320)
montage(segmentation_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'segmentations.png'));
figure(321)
montage(Arteries_fRMS_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'ArteriesfRMS.png'));
figure(3211)
montage(Vessels_velocity, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'Vessels_velocity.png'));
figure(322)
montage(ARI_velocity_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'ARIvelocity.png'));
figure(3221)
montage(VRI_velocity_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'VRIvelocity.png'));
figure(323)
montage(bvr_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'bloodVolumeRate.png'));
figure(324)
montage(histo_art_velocity_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'histogramVelocityArteries.png'));
figure(327)
montage(Stroke_total_volume, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, '_strokeAndTotalVolume.png'));
figure(328)
montage(A_sections, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'A_sections.png'));
figure(329)
montage(diasys_Artery, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'diasys_Artery.png'));
figure(330)
montage(diasys_Vein, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'diasys_Vein.png'));
figure(331)
montage(ArterialWaveformAnalysis_artery, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'ArterialWaveformAnalysis_artery.png'));

end
