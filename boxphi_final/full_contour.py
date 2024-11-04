from getdist import plots

analysis_settings = {'ignore_rows': '0.6',
                     'smooth_scale_2D': '0.4'}
g=plots.get_subplot_plotter(chain_dir=r'/data2/chains/hede_mc/final',analysis_settings=analysis_settings)
g.settings.axes_fontsize = 15
g.settings.axes_labelsize = 18
g.settings.legend_fontsize = 22
roots = []
roots.append('boxphi_f_mc_p20base')
roots.append('phi4_f_mc_p20base')
roots.append('lcdm_mc_p20base')
roots.append('boxphi_f_mc_p20des')
roots.append('phi4_f_mc_p20des')
roots.append('lcdm_mc_p20des')
params = ['ombh2', 'omch2', 'H0', 'ns', 'logA', 'tau', 'S8', 'f_hede', 'lnzc', 'xiv_phi']
g.triangle_plot(roots, params, filled=[True, True, True, False, False, False], contour_colors=['blue', 'red', 'gray', 'blue', 'red', 'gray'], contour_ls=['-', '-', '-', '--', '--', '--'], contour_lws=[1, 1, 1, 2, 2, 2], legend_labels=[r'$\mathcal{G}$EDE', 'cEDE', '$\Lambda$CDM'])
g.export()
