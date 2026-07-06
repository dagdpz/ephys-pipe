function epp_update_sorting_table(monkey,dates)
% INPUT: epp_update_sorting_table('Flaffus',[20160608,20160609])
% automatically updates '<monkey>_UltraSort.xlsx' in the sorting table path
% note that this function takes information from the "final_sorting" sheet,
% and overrides "automatic_sorting" sheet, using:
% "final_sorting" as basis and overriding with information accessable from
% Electrode_depths.m as well as sortcodes from the corresponding PLX file.
% Therefore, manually added informations (such as SNR, stability, and single ranking as well as grid hole positions) are kept
% Before the automatic information can be processed further, the user is
% asked to copy the "automatic sorting" sheet to "final sorting", to make
% sure no valuable information is lost in the automatic process

max_depth_delta_for_same_file = 50;
dag_drive=DAG_get_server_IP;
threshold_used='negthr'; % to be changed to 'both'

%main_folder             =[dag_drive filesep 'Projects' filesep 'spikesorting' filesep 'testdata' filesep 'merge_data_structure' filesep 'TDTbrain' filesep];
main_folder             =[dag_drive filesep 'Data' filesep 'TDTtanks' filesep monkey '_phys' filesep];

main_folder_content     =dir(main_folder);
main_folder_content     =main_folder_content([main_folder_content.isdir]);
main_folder_content(1:2)=[];
subfolders              ={main_folder_content.name};
if nargin>=2
    subfolders=subfolders(cellfun(@(x) ismember(str2double(x),dates),subfolders));
end


% DBpath=DAG_get_Dropbox_path;
% DBfolder=[DBpath filesep 'DAG' filesep 'phys' filesep monkey '_dpz' filesep];
DBfolder=[dag_drive 'Data' filesep 'Sorting_tables' filesep monkey filesep];
[~, sheets_available]=xlsfinfo([DBfolder  monkey(1:3) '_UltraSort.xlsx']);
if ismember('final_sorting',sheets_available)
    [~, ~, sorting_table]=xlsread([DBfolder  monkey(1:3) '_UltraSort.xlsx'],'final_sorting');
elseif ismember('automatic_sorting',sheets_available)
    [~, ~, sorting_table]=xlsread([DBfolder  monkey(1:3) '_UltraSort.xlsx'],'automatic_sorting');
else
    sorting_table={'Monkey','Session','Filenumber','Blocks','Channel','z','Unit','N_spk','Neuron_ID','Site_ID'};
end
old_table=sorting_table;
dateindex_old=DAG_find_column_index(old_table,'Session');
if ~isempty(dateindex_old)
    dates_old=[old_table{2:end,dateindex_old}];
    unique_old_dates=unique(dates_old);
else
    unique_old_dates=[];
end
sorting_table=sorting_table(1,:);


for c=1:numel(sorting_table)
    column_name = strrep(sorting_table{c},' ','_');
    column_name = strrep(column_name,'?','');
    sorting_table{c}=column_name;
    idx.(column_name)=DAG_find_column_index(sorting_table,column_name);
end
old_table(1,:)=sorting_table;

%% load electrode depths
clear Session block channels z
run([DBfolder  filesep 'Electrode_depths_' monkey(1:3)]);
Electrode_depths=struct('Session',Session,'block',block,'channels',channels,'z',z);

%% Check for apparent mistakes (nonmatching channels for same cells would lead to errors later on)
for c=1:numel(Electrode_depths)
    if numel(Electrode_depths(c).channels) ~= numel(Electrode_depths(c).z)
        disp(['Problem in ' num2str(Electrode_depths(c).Session) ', block ' num2str(Electrode_depths(c).block), ' channels and depths dont match']);
    end
end

n_row=1;
new_sessions_counter=0;
for s =1:numel(subfolders)
    date=subfolders{s};
    session=str2double(date);
    if ~ismember(session,unique_old_dates)
        new_sessions_counter=new_sessions_counter+1;
    end
    %matfiles=dir([main_folder date filesep 'spikes' filesep 'dataspikes*negthr.mat']);
    
    matfiles=dir([main_folder date filesep 'dataspikes*' threshold_used '.mat']);
    
    matfiles={matfiles.name};
    matfiles=sort(matfiles);
    unit_per_session_counter=0;
    site_per_session_counter=0;
    for f=1:numel(matfiles)
        site_per_session_counter=site_per_session_counter+1;
        matfile=matfiles{f};
        
        load([main_folder date filesep matfile]);
        %load([main_folder date filesep 'spikes' filesep matfile]);
        Channel=str2double(extractBetween(matfile,'_ch',['_' threshold_used]));
        Filenumber=str2double(extractBetween(matfile,'_rb','_ch'));
        
        
        d=ismember([Electrode_depths.Session],session);
        
        [channel_depths, channel_blocks] = get_channel_depths_and_blocks(Electrode_depths(d), Channel);
        channel_filenumbers = assign_filenumbers_from_depth(channel_depths, max_depth_delta_for_same_file);
        matching_entries = (channel_filenumbers == Filenumber);
        matched_blocks = channel_blocks(matching_entries);
        matched_depths = channel_depths(matching_entries);
        z=matched_depths(1);
        
        
        
        sitename=[monkey(1:3) '_' date  '_Site_' sprintf('%02d',site_per_session_counter)];
        units_present=unique(cluster_class(:,1));
        units_present(units_present==0)=[];
        if any(units_present)
            for sortcode=units_present'
                n_row=n_row+1;
                unit_per_session_counter=unit_per_session_counter+1;
                sorting_table{n_row,idx.z}=z;
                sorting_table{n_row,idx.Monkey}=monkey;
                sorting_table{n_row,idx.Session}=session;
                sorting_table{n_row,idx.Filenumber}=Filenumber;
                sorting_table{n_row,idx.Channel}=Channel;
                sorting_table{n_row,idx.Blocks}=serialize_blocks_for_excel(matched_blocks);
                sorting_table{n_row,idx.Unit}=sortcode;
                sorting_table{n_row,idx.N_spk}=sum(cluster_class(:,1)==sortcode);
                sorting_table{n_row,idx.Neuron_ID}=[monkey(1:3) '_' date  '_' sprintf('%03d',unit_per_session_counter)];
                sorting_table{n_row,idx.Site_ID}=sitename;
            end
        else
            sorting_table{n_row,idx.z}=z;
            sorting_table{n_row,idx.Monkey}=monkey;
            sorting_table{n_row,idx.Session}=session;
            sorting_table{n_row,idx.Filenumber}=Filenumber;
            sorting_table{n_row,idx.Channel}=Channel;
            sorting_table{n_row,idx.Blocks}=serialize_blocks_for_excel(matched_blocks);
            sorting_table{n_row,idx.Unit}=0;
            sorting_table{n_row,idx.N_spk}=0;
            sorting_table{n_row,idx.Neuron_ID}=[monkey(1:3) '_' date  '_000'];
            sorting_table{n_row,idx.Site_ID}=sitename;
        end
    end
end

%% append old information for each line in the new table
if size(old_table,1)>1
    for r=2:size(sorting_table,1)
        session_exists=[true;ismember(vertcat(old_table{2:end,idx.Session}),sorting_table{r,idx.Session})];
        site_exists= session_exists & [true; ismember(vertcat(old_table{2:end,idx.Channel}),sorting_table{r,idx.Channel}) &...
            ismember(vertcat(old_table{2:end,idx.Filenumber}),sorting_table{r,idx.Filenumber})];
        unit_exists=    site_exists & [true;ismember(vertcat(old_table{2:end,idx.Unit}),sorting_table{r,idx.Unit})];
        
        if sum(unit_exists)>1
            % this part adds unit-specific information to each line
            new_line=DAG_update_cell_table(old_table(unit_exists,:),sorting_table([1 r],:),'Session');
            sorting_table(r,:)=new_line(2,:);
        elseif sum(site_exists)>1
            % this part adds site specific information to each line
            oldrow=find(site_exists);
            oldrow=oldrow([1 2]);
            new_line=DAG_update_cell_table(old_table(oldrow,:),sorting_table([1 r],:),'Session');
            sorting_table(r,:)=new_line(2,:);
        elseif sum(session_exists)>1
            % this part adds session-specific information to each line
            oldrow=find(session_exists);
            oldrow=oldrow([1 2]);
            new_line=DAG_update_cell_table(old_table(oldrow,:),sorting_table([1 r],:),'Session');
            sorting_table(r,:)=new_line(2,:);
        end
    end
    %% remove old table entries for the sessions we are updating
    old_table([false, ismember([old_table{2:end,idx.Session}],[sorting_table{2:end,idx.Session}])],:)=[];
end
%% fill in rows corresponding to new sessions
sessions=unique([sorting_table{2:end,idx.Session}]);
complete_table=old_table;
for s=1:numel(sessions)
    session=sessions(s);
    rows_new=[false sorting_table{2:end,idx.Session}]==session;
    old_sessions=[0 complete_table{2:end,idx.Session}];
    split_after_this=find(old_sessions<session,1,'last');
    split_until_here=find(old_sessions>session,1,'first');
    if any(split_until_here)
        complete_table=[complete_table(1:split_after_this,:);...
            sorting_table(rows_new,:);...
            complete_table(split_until_here:end,:)];
    else
        complete_table=[complete_table(1:split_after_this,:);...
            sorting_table(rows_new,:)];
    end
end
% %% add empty rows at the end to overwrite excess rows in case of shorting a session
old_table_rows=size(old_table,1);
if old_table_rows>size(complete_table,1)
    complete_table=[complete_table; cell(old_table_rows-size(complete_table,1),size(complete_table,2))];
end
%[complete_mastertable]=DAG_update_cell_table(sorting_table,old_table,'Session');
xlswrite([DBfolder filesep monkey(1:3) '_UltraSort.xlsx'],complete_table,'automatic_sorting');
end

function [channel_depths, channel_blocks] = get_channel_depths_and_blocks(electrode_depths_session, channel)
channel_depths = [];
channel_blocks = [];

for i = 1:numel(electrode_depths_session)
    ch_i = electrode_depths_session(i).channels(:);
    z_i = electrode_depths_session(i).z(:);
    if numel(ch_i) ~= numel(z_i)
        continue;
    end
    keep = (ch_i == channel);
    if any(keep)
        channel_depths = [channel_depths; z_i(keep)]; %#ok<AGROW>
        channel_blocks = [channel_blocks; repmat(electrode_depths_session(i).block, sum(keep), 1)]; %#ok<AGROW>
    end
end
end

function filenumbers = assign_filenumbers_from_depth(channel_depths, max_depth_delta)
% Assign consecutive file numbers from depth progression.
% Starts at 1 and increases only when depth change exceeds max_depth_delta.

if nargin < 2 || isempty(max_depth_delta)
    error('Provide max_depth_delta as second input.');
end

if isempty(channel_depths)
    filenumbers = [];
    return;
end

is_row = isrow(channel_depths);
depths = channel_depths(:);
filenumbers = zeros(size(depths));
filenumbers(1) = 1;

for i = 2:numel(depths)
    prev_depth = depths(i-1);
    curr_depth = depths(i);
    
    % If one of the two depths is invalid, start a new file number.
    if ~isfinite(prev_depth) || ~isfinite(curr_depth)
        filenumbers(i) = filenumbers(i-1) + 1;
    elseif abs(curr_depth - prev_depth) <= max_depth_delta
        filenumbers(i) = filenumbers(i-1);
    else
        filenumbers(i) = filenumbers(i-1) + 1;
    end
end

if is_row
    filenumbers = filenumbers.';
end
end

function blocks_out = serialize_blocks_for_excel(blocks_in)
% Convert block values to Excel-safe text representation.
% Non-empty values are always stored as comma-separated text (e.g. '12,14').

if isempty(blocks_in)
    blocks_out = '';
    return;
end

blocks_vec = unique(blocks_in(:).');
blocks_vec = blocks_vec(isfinite(blocks_vec));
if isempty(blocks_vec)
    blocks_out = '||';
else
    blocks_out = sprintf('%g|', blocks_vec);
    blocks_out = ['|' blocks_out];
end
end

