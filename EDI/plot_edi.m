function fig = plot_edi(time_vec, edi, edi_stats, cat_edi, cat_edi_stats, cat_names, categories_vec)
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

plot_cdi = false; plot_edi_per_cat = false; show_stats = false; mask = [];
fig_sz = [0.15 0.15 0.25 0.3]; marg_w = [0.085, 0.05]; marg_h = [0.14 0.08];
lims.y = get_ylims(edi, 20);
make_nice = []; make_nice.keep_ticks = true; make_nice.line_col = linspecer(1); make_nice.stim_len = stim_len;
    
if isfield(cat_edi,'cdi')
    plot_cdi = true;
    fig_sz(4) = 0.4; marg_h(1) = 0.11;
    ylims2 = get_ylims(cat_edi.cdi, 20);
    lims.y(1) = min([lims.y(1) ylims2(1)]);
    lims.y(2) = max([lims.y(2) ylims2(2)]);
    make_nice.keep_ticks = false; 
end
if isfield(cat_edi,'edi_per_categ')
    plot_edi_per_cat = true;
    n_cat = size(cat_edi.edi_per_categ, 2);
    fig_sz(3) = 0.4; marg_w = [0.05 0.5];
end
if isfield(edi_stats,'p_values')
    show_stats = true;
end

fig = figure('Units','Normalized','Position',fig_sz);
ha = tight_subplot(1+double(plot_cdi), 1, 0.06, marg_h, marg_w);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
if show_stats; mask = edi_stats.mask; end
tit_info = []; tit_info.txt = 'EDI (Exemplar discriminability index)'; tit_info.change = false;
single_edi_plot(ha(1), edi, time_vec, mask, tit_info, lims, make_nice)
if plot_cdi
    if show_stats; mask = cat_edi_stats.cdi.mask; end
    tit_info.txt = 'CDI (Category discriminability index)';
    make_nice.keep_ticks = true;
    single_edi_plot(ha(2), cat_edi.cdi, time_vec, mask, tit_info, lims, make_nice)
end

if plot_edi_per_cat
    ha2 = tight_subplot(n_cat, 1, 0.05, [0.11 0.08], [0.55, 0.03]);
    % Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right] 
    sup_ax = tight_subplot(1, 1, 0.05, marg_h, [0.55, 0.03]);
    title(sup_ax, 'Single category EDIs', 'Visible', 'on'); set(sup_ax, 'Visible', 'off')
    lims.y = get_ylims(cat_edi.edi_per_categ, 20); make_nice.keep_ticks = false;
    linecols = linspecer(n_cat); tit_info.change = true;
    for categ = 1:n_cat
        make_nice.line_col = linecols(categ,:);
        if show_stats; mask = cat_edi_stats.edi_per_categ.mask(:, categ); end
        if categ == n_cat; make_nice.keep_ticks = true; end
        tit_info.txt = sprintf('%s (n_{stim}=%d)',cat_names{categ},sum(categories_vec==categ));
        single_edi_plot(ha2(categ), cat_edi.edi_per_categ(:, categ), time_vec, mask, tit_info, lims, make_nice)
    end
end

end

function single_edi_plot(ax_handle, dat, time_vec, mask, tit_info, lims, make_nice)
axes(ax_handle);
fig_col = get(gcf, 'Color');
mask_dat = nan(size(time_vec)); mask_dat(boolean(mask)) = 0;

plot(time_vec, dat, 'LineWidth', 2, 'Color', make_nice.line_col); hold on
plot([0 0], lims.y, 'k:', 'LineWidth', 0.5);
plot(lims.x, [0 0], 'k:', 'LineWidth', 0.5);
plot([make_nice.stim_len make_nice.stim_len], lims.y, 'k:', 'LineWidth', 0.5);
plot(time_vec, mask_dat, 'r-', 'LineWidth',2.5);

xlim(lims.x); ylim(lims.y)
tit_handle = title(tit_info.txt, 'FontSize', 11); set(tit_handle,'Units','Normalized');
if tit_info.change
    set(tit_handle,'FontSize',10,'Position',[1,1.15,0],'FontWeight','Normal','HorizontalAlignment','right','VerticalAlignment','top');
end
if make_nice.keep_ticks
    xlabel('Time (ms)');
else
    set(gca,'YTickLabels','','XTickLabels','');
end
set(gca,'FontSize',10,'Color',fig_col); box off

end

function ylims = get_ylims(dat, rounder)
ylims = [floor(min(dat(:))*rounder)/rounder, ceil(max(dat(:))*rounder)/rounder];
end