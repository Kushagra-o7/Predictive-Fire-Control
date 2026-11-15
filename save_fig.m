function save_fig(fig, base_path)
%SAVE_FIG  Save a figure as PNG (for the paper) and .fig (editable in MATLAB).
folder = fileparts(base_path);
if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
try
    exportgraphics(fig, [base_path '.png'], 'Resolution', 200);
catch
    saveas(fig, [base_path '.png']);          % older MATLAB / Octave
end
try
    savefig(fig, [base_path '.fig']);
catch
end
end
