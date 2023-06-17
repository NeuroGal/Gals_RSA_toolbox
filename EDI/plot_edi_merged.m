function fig = plot_edi_merged(time_vec, edi, edi_stats, cat_edi, cat_edi_stats, cat_names, categories_vec)
% edi, edi_stats, cat_edi, cat_edi_stats are all outputs of edi_calc
%
% Gal Vishne, 2022, gal.vishne@gmail.com
%
% Uses tight_subplot & linspecer:
%   https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
%   https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap

stim_len = 0; % change here if you want a vertical line somewhere
lims = [];
lims.x = [time_vec(1) time_vec(end)];
show_stats = isfield(edi_stats,'p_values');
mask = [];
rounder = 10;

plot_cdi = isfield(cat_edi,'cdi'); plot_edi_per_cat = isfield(cat_edi,'edi_per_categ'); 
nrow =  1+double(plot_cdi)+double(plot_edi_per_cat);

fig_sz = [0.1 0.1 0.25 0.1+nrow*0.15]; marg_w = [0.085, 0.05]; marg_h = [0.1 0.08];
make_nice = []; make_nice.line_col = 0.3*ones(3,1); make_nice.stim_len = stim_len;
make_nice.mask_col = 'r'; make_nice.mask_offset = 0.03; make_nice.line_width = 2;
make_nice_defaults = make_nice;

fig = figure('Units','Normalized','Position',fig_sz);
ha = tight_subplot(nrow, 1, 0.16, marg_h, marg_w);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]

if show_stats; mask = edi_stats.mask; end
lims.y = get_ylims(edi, rounder);
tit = 'EDI (Exemplar discriminability index)';
single_edi_plot(ha(1), edi, time_vec, mask, tit, lims, make_nice_defaults);
if plot_edi_per_cat
    tit = 'Single category EDIs'; n_cat = size(cat_edi.edi_per_categ, 2);
    lims.y = get_ylims(cat_edi.edi_per_categ, rounder); % lims.y(2) = lims.y(2) + 0.05;
    linecols = linspecer(n_cat+1); linecols = linecols([1 3:(n_cat+1)],:);
    legtxt = cell(n_cat,1); plts = [];
    for categ = 1:n_cat
        make_nice.line_col = [linecols(categ,:) 0.7]; %alpha is 4th arg
        make_nice.line_width = 1.5;
        if show_stats
            mask = cat_edi_stats.edi_per_categ.mask(:, categ); 
            make_nice.mask_offset = 0.0125*categ;
            make_nice.mask_col = make_nice.line_col;
        end
        if all(isnan(cat_edi.edi_per_categ(:, categ)))
        else
            legtxt{categ} = sprintf('%s (%d images)',cat_names{categ},sum(categories_vec==categ));
            plts(categ) = single_edi_plot(ha(2), cat_edi.edi_per_categ(:, categ), time_vec, mask, tit, lims, make_nice);
        end
    end
    make_nice = make_nice_defaults; make_nice.mask_offset = 0.0125*(n_cat+1); make_nice.line_width = 2.5;
    legtxt{n_cat+1} = 'Mean single categories';
    plts(n_cat+1) = single_edi_plot(ha(2), cat_edi.edi_per_categ_mean, time_vec, cat_edi_stats.edi_per_categ_mean.mask, tit, lims, make_nice);
    leg = legend(ha(2), plts(~cellfun(@isempty,legtxt)), legtxt(~cellfun(@isempty,legtxt)), 'Box', 'off');
    leg.Position(1:2) = [0.625 0.56];
end
if plot_cdi
    lims.y = get_ylims(cat_edi.cdi, rounder);
    if show_stats; mask = cat_edi_stats.cdi.mask; end
    tit = 'CDI (Category discriminability index)';
    single_edi_plot(ha(3), cat_edi.cdi, time_vec, mask, tit, lims, make_nice_defaults);
end

end

function plt_handle = single_edi_plot(ax_handle, dat, time_vec, mask, tit, lims, make_nice)
axes(ax_handle);
fig_col = get(gcf, 'Color');
mask_dat = nan(size(time_vec)); mask_dat(boolean(mask)) = lims.y(1)+make_nice.mask_offset;

plt_handle = plot(time_vec, dat, 'LineWidth', make_nice.line_width, 'Color', make_nice.line_col); hold on
plot([0 0], lims.y, 'k:', 'LineWidth', 0.5);
plot(lims.x, [0 0], 'k:', 'LineWidth', 0.5);
plot([make_nice.stim_len make_nice.stim_len], lims.y, 'k:', 'LineWidth', 0.5);
plot(time_vec, mask_dat, 'Color', make_nice.mask_col, 'LineWidth', make_nice.line_width+0.5);

xlim(lims.x); ylim(lims.y)
tit_handle = title(tit, 'FontSize', 11,'Units','Normalized');
tit_handle.HorizontalAlignment = 'left'; tit_handle.Position(1) = 0.05;
xlabel('Time (ms)');
set(gca,'FontSize',10,'Color',fig_col); box off
end

function ylims = get_ylims(dat, rounder)
ylims = [floor(min(dat(:))*rounder)/rounder, ceil(max(dat(:))*rounder)/rounder];
end