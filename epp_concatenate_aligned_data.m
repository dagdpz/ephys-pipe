function combined = epp_concatenate_aligned_data(combined, chunk)
combined.trial = [combined.trial; chunk.trial(:)];                                          % #ok<AGROW>
combined.state.Value = [combined.state.Value; chunk.state.Value(:)];                        % #ok<AGROW>
combined.state.timestamps = [combined.state.timestamps; chunk.state.timestamps(:)];         % #ok<AGROW>
combined.sample.x_hnd = [combined.sample.x_hnd; chunk.sample.x_hnd(:)];                     % #ok<AGROW>
combined.sample.y_hnd = [combined.sample.y_hnd; chunk.sample.y_hnd(:)];                     % #ok<AGROW>
combined.sample.x_eye = [combined.sample.x_eye; chunk.sample.x_eye(:)];                     % #ok<AGROW>
combined.sample.y_eye = [combined.sample.y_eye; chunk.sample.y_eye(:)];                     % #ok<AGROW>
combined.sample.timestamps = [combined.sample.timestamps; chunk.sample.timestamps(:)];      % #ok<AGROW>
end
