function epp_update_sorting_table(monkey,dates,cfg)
%EPP_UPDATE_SORTING_TABLE  Rebuild the sorting sheet from spikes.
%
%   epp_update_sorting_table(monkey)
%   epp_update_sorting_table(monkey, dates)
%   epp_update_sorting_table(monkey, dates, cfg)
%
%   Inputs:
%     monkey - char/string, at least 3 characters (e.g. 'Flaffus').
%     dates  - optional numeric vector of YYYYMMDD session dates. Omit to scan
%              every YYYYMMDD subfolder.
%
%   Reads (not via cfg; paths come from DAG_get_server_IP):
%     <drive>\Data\TDTtanks\<monkey>_phys\<YYYYMMDD>\dataspikes_*<spike_threshold>.mat
%         cluster_class col 1 = sort code, col 2 = spike time in ms.
%     <drive>Data\Sorting_tables\<monkey>\Electrode_depths_<monkey(1:3)>.m
%         script that must assign Session, block, channels, z.
%     <drive>Data\Sorting_tables\<monkey>\<monkey(1:3)>_UltraSort.xlsx
%         sheet cfg.sorting_sheet only (default automatic_sorting). Missing file starts from a header row.
%         Missing sheet is an error (a backup is still written).
%
%   Outputs:
%     none. Before writing, copies the workbook to
%     <Mon>_UltraSort_backup_yyyymmdd_HHMMSS.xlsx in the same folder.
%     Then replaces the workbook with a single sheet, cfg.sorting_sheet.
%     Columns filled here: Monkey, Session, Filenumber, Blocks, Channel, z,
%     Unit, N_spk, Neuron_ID, Site_ID. Other columns are copied from the
%     previous automatic_sorting sheet when session/channel/file/unit match.
%     Blocks are pipe text ('|2|5|'). Unit 0 marks a channel with no sort codes.
%
% INPUT: epp_update_sorting_table('Flaffus',[20160608,20160609])
% Rebuilds automatic_sorting from Electrode_depths and spike sort codes.
% Manual columns (SNR, stability, Dataset, Hemisphere, Perturbation, ...)
% are kept when the previous automatic_sorting row matches.

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(genpath(fullfile(repo_root, 'external')));
max_depth_delta_for_same_file = 50;
dag_drive=DAG_get_server_IP;
if nargin < 3 || isempty(cfg)
    cfg = epp_general_settings('_', struct('project_version', ''));
end
threshold_used = cfg.spike_threshold;
sheetname = cfg.sorting_sheet;

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
xlsx_path = [DBfolder monkey(1:3) '_UltraSort.xlsx'];
if isfile(xlsx_path)
    backup_sorting_workbook(xlsx_path);
    [~, sheets_available] = xlsfinfo(xlsx_path);
    if ~ismember(sheetname, sheets_available)
        error('epp_update_sorting_table:MissingSheet', ...
            '%s has no %s sheet. A backup was written next to the workbook. Copy your table onto a sheet named %s and rerun.', ...
            xlsx_path, sheetname, sheetname);
    end
    [~, ~, sorting_table] = xlsread(xlsx_path, sheetname);
else
    sorting_table = empty_sorting_header();
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
        
        sss=load([main_folder date filesep matfile]);
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
            % Empty channel is its own site row. Unit 0, Neuron_ID *_000.
            n_row=n_row+1;
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
%[complete_mastertable]=DAG_update_cell_table(sorting_table,old_table,'Session');
write_automatic_sorting_workbook(xlsx_path, complete_table, sheetname);
end

function sorting_table = empty_sorting_header()
sorting_table = {'Monkey','Session','Filenumber','Blocks','Channel','z','Unit','N_spk','Neuron_ID','Site_ID'};
end

function backup_sorting_workbook(xlsx_path)
% Full copy of the workbook before automatic_sorting is replaced.
[folder, name, ext] = fileparts(xlsx_path);
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
backup_path = fullfile(folder, sprintf('%s_backup_%s%s', name, stamp, ext));
suffix = 1;
while isfile(backup_path)
    suffix = suffix + 1;
    backup_path = fullfile(folder, sprintf('%s_backup_%s_%d%s', name, stamp, suffix, ext));
end
copyfile(xlsx_path, backup_path);
end

function write_automatic_sorting_workbook(xlsx_path, complete_table, sheet_name)
% Write a new workbook whose only sheet is sheet_name, then replace
% the live file. The pre-update backup is already on disk.
[folder, name, ext] = fileparts(xlsx_path);
temp_path = fullfile(folder, [name '_writing' ext]);
if isfile(temp_path)
    delete(temp_path);
end
xlswrite(temp_path, complete_table, sheet_name);
keep_only_sheet(temp_path, sheet_name);
if isfile(xlsx_path)
    delete(xlsx_path);
end
movefile(temp_path, xlsx_path);
end

function keep_only_sheet(xlsx_path, sheet_name)
% xlswrite on a new file leaves the default Sheet1..Sheet3 behind.
excel = actxserver('Excel.Application');
excel.DisplayAlerts = false;
excel.Visible = false;
wb = excel.Workbooks.Open(xlsx_path);
cleanup_obj = onCleanup(@() release_excel(excel, wb)); %#ok<NASGU>
sheets = wb.Sheets;
has_target = false;
for i = 1:sheets.Count
    if strcmp(sheets.Item(i).Name, sheet_name)
        has_target = true;
        break;
    end
end
if ~has_target
    sheets.Item(1).Name = sheet_name;
end
for i = sheets.Count:-1:1
    if ~strcmp(sheets.Item(i).Name, sheet_name)
        sheets.Item(i).Delete;
    end
end
wb.Save;
end

function release_excel(excel, wb)
try
    if ~isempty(wb)
        wb.Close(false);
    end
catch
end
excel.Quit;
delete(excel);
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

