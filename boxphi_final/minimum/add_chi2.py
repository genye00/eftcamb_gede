import numpy as np
import os

for file in os.listdir('./'):
    if '.txt' in file:
        fin = open(file, 'r')
        nms = fin.readline().replace('#','').split()
        fin.close()
        dat = np.loadtxt(file)
        dat[nms.index('chi2__CMB')]=dat[nms.index('chi2__planck_2020_hillipop.TTTEEE')]+dat[nms.index('chi2__planck_2020_lollipop.lowlE')]+dat[nms.index('chi2__planck_2018_lowl.TT')]
        with open('./tmp/'+file, 'w') as fou:
            fou.write('#\t' + '\t'.join(nms) + '\n')
            for i,nm in enumerate(nms):
                if nm=='weight':
                    fou.write('\t%d'%dat[i])
                else:
                    fou.write('\t%.8e'%dat[i])
            fou.write('\n')


